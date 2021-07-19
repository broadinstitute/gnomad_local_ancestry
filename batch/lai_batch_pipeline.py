# noqa: D100
import logging
import hailtop.batch as hb

from gnomad.utils.slack import slack_notifications
from batch.batch_utils import (
    init_arg_parser,
    run_batch,
)

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

GCLOUD_USER_ACCOUNT = "mwilson@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://seqr-mwilson-misc/"


def check_args(parser, args):
    """
    Check passed args to ensure pipeline can run properly.

    :param parser: arg parser.
    :param args: Arg's from argparser.
    """
    if not (args.run_eagle or args.run_rfmix or args.run_tractor):
        parser.error(
            "Need to specify at least one tool to run (--run-eagle, --run-rfmix, and/or --run-tractor)"
        )
    if args.run_eagle:
        if not (args.sample_vcf or args.ref_vcf):
            parser.error(
                "Need to specify either sample and/or reference vcfs (--sample-vcf or --ref-vcf)"
            )
    if args.run_rfmix:
        if not ((args.run_eagle and args.sample_vcf) or args.phased_sample_vcf):
            parser.error(
                "Need to specify either sample vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --sample-vcf or --phased-sample-vcf)"
            )
        if not ((args.run_eagle and args.ref_vcf) or args.phased_ref_vcf):
            parser.error(
                "Need to specify either reference vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --reference-vcf or --phased-reference-vcf)"
            )
        if not args.genetic_map:
            parser.error(
                "Need to specify genetic recombination map, --genetic-map, for RFMix to run."
            )
        if not args.pop_sample_map:
            parser.error(
                "Need to specify sample to population mapping file, --pop-sample-map"
            )

    if args.run_tractor:
        if not ((args.run_eagle and args.sample_vcf) or args.phased_sample_vcf):
            parser.error(
                "Need to specify either sample vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --sample-vcf or --phased-sample-vcf)"
            )
        if not args.run_rfmix and not args.msp_file:
            parser.error(
                "Need to run RFMix to generate MSP file or pass the MSP tsv file to the script, --msp-file."
            )
        if not args.ancs:
            parser.error(
                "Need to specify either number of continental ancestries within RFMix and phased sample VCF."
            )


def eagle(
    batch,
    vcf,
    contig,
    mem="standard",
    storage="50G",
    threads=8,
    image="gcr.io/broad-mpg-gnomad/lai_phasing:latest",
) -> hb.Batch.new_job:
    """
    Run the phasing tool Eagle on passed VCF.

    :param batch: Hail batch
    :param vcf: VCF to phase
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "standard".
    :param storage: Hail batch job storage, deaults to "50G".
    :param threads: The number of threads, should match the number of CPUs requested, defaults to 8.
    :param image: Docker image for eagle job, defaults to "gcr.io/broad-mpg-gnomad/lai_phasing:latest".
    :return: Batch job
    """
    e = batch.new_job(name=f"Eagle - chr{contig}")
    e.cpu(threads)
    e.memory(mem)
    e.storage(storage)
    e.image(image)
    e.declare_resource_group(ofile={"vcf": "{root}.vcf.gz"})

    cmd = f"""
    Eagle_v2.4.1/eagle \
        --geneticMapFile Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
        --numThreads {threads} \
        --outPrefix {e.ofile} \
        --vcfOutFormat z \
        --vcf {vcf}
    """

    e.command(cmd)
    return e


def rfmix(
    batch,
    sample_pvcf,
    ref_pvcf,
    contig,
    sample_map,
    rf_genetic_map,
    mem="highmem",
    storage="50G",
    threads=8,
    image="gcr.io/broad-mpg-gnomad/lai_rfmix:latest",
) -> hb.Batch.new_job:
    """
    Run RFMix2 on phased VCF.

    :param batch: Hail batch object
    :param sample_pvcf: Phased sample VCF from phasing tool like Eagle or SHAPEIT.
    :param ref_pvcf: Phased reference sample VCF from phasing tool like Eagle or ShapeIt.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param sample_map: TSV file containing a mapping from sample IDs to ancestral populations, i.e. NA12878    EUR.
    :param rf_genetic_map: HapMap genetic map from SNP base pair positions to genetic coordinates in cM.
    :param mem: Hail batch job memory, defaults to "standard".
    :param storage: Hail batch job storage, deaults to "50G".
    :param threads: The number of threads, should match the number of CPUs requested, defaults to 8.
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_tractor:latest"
    :return: Hail batch job
    """
    r = batch.new_job(name=f"RFMix - chr{contig}")
    r.memory(mem)
    r.storage(storage)
    r.cpu(threads)
    r.image(image)
    r.declare_resource_group(
        ofile={"msp.tsv": "{root}.msp.tsv", "fb.tsv": "{root}.fb.tsv"}
    )
    cmd = f"""
    ./rfmix/rfmix \
        -f {sample_pvcf} \
        -r {ref_pvcf} \
        --chromosome=chr{contig} \
        -m {sample_map} \
        -g {rf_genetic_map} \
        -n 5 \
        -e 1 \
        --reanalyze-reference \
        -o {r.ofile}
    """
    r.command(cmd)
    return r


def tractor(
    batch,
    msp,
    vcf,
    ancs,
    zipped,
    contig,
    mem="standard",
    storage="50G",
    image="gcr.io/broad-mpg-gnomad/lai_tractor:latest",
) -> hb.Batch.new_job:
    """
    Run Tractor's ExtractTract.py script.

    :param batch: Hail batch object
    :param msp: MSP tsv file from LAI tool like RFMix2 or XGMix.
    :param vcf: Phased sample VCF from phasing tool like Eagle or SHAPEIT.
    :param ancs: Number of ancestral population within the MSP file.
    :param zipped: Whether the input VCF file is zipped or not, i.e. ends in vcf.gz.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "standard".
    :param storage: Hail batch job storage, deaults to "50G".
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_tractor:latest"
    :return: Hail Batch job
    """
    t = batch.new_job(name=f"Tractor - chr{contig}")
    t.storage(storage)
    t.image(image)
    t.memory(mem)

    rg_def = {}
    for i in range(ancs):
        rg_def[f"vcf{i}"] = f"{{root}}.anc{i}.vcf"
        rg_def[f"dos{i}.txt"] = f"{{root}}.anc{i}.dosage.txt"
        rg_def[f"ancdos{i}.txt"] = f"{{root}}.anc{i}.hapcount.txt"

    t.declare_resource_group(ofile=rg_def)
    zipped = "--zipped" if zipped else ""
    cmd = f"""
        python3 Tractor/ExtractTracts.py --msp {msp} --vcf {vcf} --num-ancs={ancs} {zipped} --output-path={t.ofile}
        """
    t.command(cmd)
    t.command(f"ls -l {t.ofile}*")
    return t


def run_lai(args):
    """Run batch LAI pipeline."""
    contig = args.contig
    logger.info(f"Running gnomAD LAI on chr{contig}")
    with run_batch(args, f"LAI - chr{contig}") as b:
        output_path = args.output_bucket

        if args.run_eagle:
            vcf = b.read_input(args.sample_vcf)

            if args.ref_vcf:
                ref_vcf = b.read_input(args.ref_vcf)
                ref_e = eagle(b, ref_vcf, contig)
                b.write_output(
                    ref_e.ofile, dest=f"{output_path}eagle/chr{contig}_reference"
                )

            e = eagle(b, vcf, contig)
            b.write_output(
                e.ofile, dest=f"{output_path}eagle/output/phased_chr{contig}_amr"
            )

        if args.run_rfmix:
            sample_map = b.read_input(args.pop_sample_map)
            genetic_map = b.read_input(args.genetic_map)
            if args.phased_ref_vcf:
                phased_ref_vcf = b.read_input(args.phased_ref_vcf)
            else:
                phased_ref_vcf = ref_e.ofile
            if args.phased_sample_vcf:
                phased_sample_vcf = b.read_input(args.phased_sample_vcf)
            else:
                phased_sample_vcf = e.ofile.vcf

            r = rfmix(
                b, phased_sample_vcf, phased_ref_vcf, contig, sample_map, genetic_map
            )
            b.write_output(r.ofile, dest=f"{output_path}rfmix/output/lai_chr20_amr")

        if args.run_tractor:
            if args.phased_sample_vcf:
                phased_sample_vcf = b.read_input_group(
                    **{"vcf.gz": args.phased_sample_vcf}
                )
            else:
                phased_sample_vcf = e.ofile.vcf
            if args.msp_file:
                msp_file = b.read_input_group(**{"msp.tsv": args.msp_file})
            else:
                msp_file = r.ofile
            t = tractor(
                b, msp_file, phased_sample_vcf, args.ancs, zipped=True, contig=contig
            )
            b.write_output(t.ofile, dest=f"{output_path}tractor/output/chr20_amr2_11s")


if __name__ == "__main__":
    p = init_arg_parser(
        default_cpu=8,
        default_billing_project="broad-mpg-gnomad",
        default_temp_bucket="gnomad-batch",
    )
    p.add_argument(
        "--sample-vcf",
        required=False,
        help="Google bucket path to sample VCF to phase.",
    )
    p.add_argument(
        "--ref-vcf",
        required=False,
        help="Google bucket path reference VCF to phase if separate.",
    )
    p.add_argument(
        "--output-bucket",
        required=True,
        help="Google bucket path for results. Each tool will create a subfolder here.",
    )
    p.add_argument(
        "--run-rfmix",
        required=False,
        action="store_true",
        help="Whether to run RFMix2.",
    )
    p.add_argument(
        "--genetic-map",
        required=False,
        help="Genetic map for required for RFMix2.",
        default="gs://gnomad-batch/mwilson/lai/inputs/rfmix/genetic_map_hg38.txt",
    )
    p.add_argument(
        "--pop-sample-map", required=False, help="Sample population mapping for RFMix2."
    )
    p.add_argument("--contig", required=True, help="Chromosomes to run LAI on.")
    p.add_argument(
        "--slack-channel",
        required=False,
        help="Slack channel to send job status to, needs @ for DM.",
    )
    p.add_argument(
        "--phased-ref-vcf",
        required=False,
        help="Phased reference VCF, if supplied, will not re-run phasing.",
    )
    p.add_argument(
        "--run-eagle",
        reuqired=False,
        action="store_true",
        help="Whether to run eagle to phase samples.",
    )
    p.add_argument("--phased-sample-vcf", help="VCF of phased samples.")
    p.add_argument(
        "--run-tractor",
        required=False,
        action="store_true",
        help="Run Tractor's ExtractTracts.py script.",
    )
    p.add_argument(
        "--msp-file", required=False, help="Output from LAI program like RFMix2."
    )
    p.add_argument(
        "--ancs",
        required=False,
        help="Number of ancestries within the reference panel. Used to extract ancestry tracts from phased VCF in Tractor.",
        default=3,
        type=int,
    )
    args = p.parse_args()
    check_args(args)

    if args.slack_channel:
        from slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            run_lai(args)
    else:
        run_lai(args)

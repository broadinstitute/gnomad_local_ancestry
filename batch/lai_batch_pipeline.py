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


def check_args(parser, args):
    """
    Check passed args to ensure pipeline can run properly.

    :param parser: Arg parser.
    :param args: Args from argparser.
    """
    if not (
        args.run_eagle
        or args.run_rfmix
        or args.run_xgmix
        or args.run_tractor
        or args.make_lai_vcf
    ):
        parser.error(
            "Need to specify at least one tool to run (--run-eagle, --run-rfmix, --run-xgmix, --run-tractor, and/or --make-lai-vcf)."
        )
    if args.run_eagle:
        if not (args.sample_vcf or args.ref_vcf):
            parser.error(
                "Need to specify either sample and/or reference vcfs (--sample-vcf or --ref-vcf)."
            )
    if args.run_rfmix and args.run_xgmix:
        parser.error(
            "Can only specify one LAI tool, either RFMix or XGMix (--run-rfmix or --run-xgmix)."
        )
    if args.run_rfmix or args.run_xgmix:
        if not ((args.run_eagle and args.sample_vcf) or args.phased_sample_vcf):
            parser.error(
                "Need to specify either sample vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --sample-vcf or --phased-sample-vcf)."
            )
        if not ((args.run_eagle and args.ref_vcf) or args.phased_ref_vcf):
            parser.error(
                "Need to specify either reference vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --reference-vcf or --phased-reference-vcf)."
            )
        if not args.genetic_map:
            parser.error(
                "Need to specify genetic recombination map, --genetic-map, for RFMix to run."
            )
        if not args.pop_sample_map:
            parser.error(
                "Need to specify sample to population mapping file, --pop-sample-map."
            )

    if args.run_tractor:
        if not ((args.run_eagle and args.sample_vcf) or args.phased_sample_vcf):
            parser.error(
                "Need to specify either sample vcf for eagle to run or pass a phased sample vcf for RFMix to run (--run-eagle and --sample-vcf or --phased-sample-vcf)."
            )
        if not (args.run_rfmix or args.run_xgmix) and not args.msp_file:
            parser.error(
                "Need to run RFMix to generate MSP file or pass the MSP tsv file to the script, --msp-file."
            )
        if not args.ancs:
            parser.error(
                "Need to specify either number of continental ancestries within RFMix and phased sample VCF."
            )


def eagle(
    batch: hb.Batch,
    vcf: str,
    contig: str,
    mem: str = "standard",
    storage: str = "50G",
    threads: int = 8,
    image: str = "gcr.io/broad-mpg-gnomad/lai_phasing:latest",
) -> hb.Batch.new_job:
    """
    Run the phasing tool Eagle on passed VCF.

    :param batch: Hail batch object.
    :param vcf: VCF to phase.
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
    batch: hb.Batch,
    sample_pvcf: str,
    ref_pvcf: str,
    contig: str,
    sample_map: str,
    rf_genetic_map: str,
    mem: str = "highmem",
    storage: str = "50G",
    threads: int = 8,
    image: str = "gcr.io/broad-mpg-gnomad/lai_rfmix:latest",
) -> hb.Batch.new_job:
    """
    Run RFMix2 on phased VCF.

    :param batch: Hail batch object.
    :param sample_pvcf: Phased sample VCF from phasing tool like Eagle or SHAPEIT.
    :param ref_pvcf: Phased reference sample VCF from phasing tool like Eagle or SHAPEIT.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param sample_map: TSV file containing a mapping from sample IDs to ancestral populations, i.e. NA12878    EUR.
    :param rf_genetic_map: HapMap genetic map from SNP base pair positions to genetic coordinates in cM.
    :param mem: Hail batch job memory, defaults to "standard".
    :param storage: Hail batch job storage, deaults to "50G".
    :param threads: The number of threads, should match the number of CPUs requested, defaults to 8.
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_rfmix:latest".
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
    ./rfmix \
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


def xgmix(
    batch: hb.Batch,
    sample_pvcf: str,
    xg_genetic_map: str,
    contig: str,
    ref_pvcf: str,
    sample_map: str,
    mem: str = "30Gi",
    storage: str = "100G",
    threads: int = 8,
    image: str = "gcr.io/broad-mpg-gnomad/lai_xgmix:latest",
) -> hb.Batch.new_job:
    """
    Run XGMix on phased VCF.

    :param batch: Hail batch object.
    :param sample_pvcf: Phased sample VCF from phasing tool like Eagle or SHAPEIT.
    :param xg_genetic_map: HapMap genetic map from SNP base pair positions to genetic coordinates in cM.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param ref_pvcf: Phased reference sample VCF from phasing tool like Eagle or SHAPEIT.
    :param sample_map: TSV file containing a mapping from sample IDs to ancestral populations, i.e. NA12878    EUR.
    :param mem: Hail batch job memory, defaults to "standard".
    :param storage: Hail batch job storage, deaults to "50G".
    :param threads: The number of threads, should match the number of CPUs requested, defaults to 8.
    :param image: XGMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_xgmix:latest".
    :return: Hail batch job
    """
    x = batch.new_job(name=f"XGMix - chr{contig}")
    x.memory(mem)
    x.storage(storage)
    x.cpu(threads)
    x.image(image)
    x.declare_resource_group(
        ofile={"msp.tsv": "{root}.msp.tsv", "fb.tsv": "{root}.fb.tsv"}
    )
    cmd = f"""
    python3 XGMIX.py {sample_pvcf} {xg_genetic_map} /io/tmp/xgmix/output chr{contig} False {ref_pvcf} {sample_map}
    ln -s /io/tmp/xgmix/output/output.msp.tsv {x.ofile['msp.tsv']}
    ln -s /io/tmp/xgmix/output/output.fb.tsv {x.ofile['fb.tsv']}
    """
    x.command(cmd)
    return x


def tractor(
    batch: hb.Batch,
    msp: str,
    vcf: str,
    ancs: int,
    zipped: bool,
    zip_output: bool,
    contig: str,
    mem: str = "highmem",
    storage: str = "200G",
    image: str = "gcr.io/broad-mpg-gnomad/lai_tractor:latest",
) -> hb.Batch.new_job:
    """
    Run Tractor's ExtractTract.py script.

    :param batch: Hail batch object.
    :param msp: MSP tsv file from LAI tool like RFMix2 or XGMix.
    :param vcf: Phased sample VCF from phasing tool like Eagle or SHAPEIT.
    :param ancs: Number of ancestral population within the MSP file.
    :param zipped: Whether the input VCF file is zipped or not, i.e. ends in vcf.gz.
    :param zip_output: Whether to zip the tool's output files.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "200G".
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_tractor:latest".
    :return: Hail Batch job
    """
    t = batch.new_job(name=f"Tractor - chr{contig}")
    t.storage(storage)
    t.image(image)
    t.memory(mem)

    rg_def = {}
    file_extension = {".gz" if zip_output else ""}
    for i in range(ancs):
        rg_def[f"vcf{i}{file_extension}"] = f"{{root}}.anc{i}.vcf{file_extension}"
        rg_def[
            f"dos{i}.txt{file_extension}"
        ] = f"{{root}}.anc{i}.dosage.txt{file_extension}"
        rg_def[
            f"ancdos{i}.txt{file_extension}"
        ] = f"{{root}}.anc{i}.hapcount.txt{file_extension}"

    t.declare_resource_group(ofile=rg_def)
    zipped = "--zipped" if zipped else ""
    zip_output = "--zip-output" if zip_output else ""
    cmd = f"""
        python3 ExtractTracts.py --msp {msp} --vcf {vcf} --num-ancs={ancs} {zipped} {zip_output} --output-path={t.ofile}
        """
    t.command(cmd)
    return t


def generate_lai_vcf(
    batch: hb.Batch,
    msp: str,
    tractor_output: str,
    zipped: bool,
    contig: str,
    mem: str = "highmem",
    storage: str = "200G",
    image: str = "gcr.io/broad-mpg-gnomad/lai_vcf:latest",
) -> hb.Batch.new_job:
    """
    Run generate_output_vcf.py script.

    :param batch: Hail batch object.
    :param msp: MSP tsv file from LAI tool like RFMix2 or XGMix.
    :param tractor_output: Path to Tractor's output files.
    :param ancs: Number of ancestral populations within the MSP file.
    :param zipped: Whether the input VCF file is zipped or not, i.e. ends in vcf.gz.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "200G".
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_tractor:latest".
    :return: Hail Batch job
    """
    v = batch.new_job(name=f"Generate final VCF - chr{contig}")
    v.storage(storage)
    v.image(image)
    v.memory(mem)
    v.declare_resource_group(ofile={"vcf.bgz": "{root}_lai_annotated.vcf.bgz"})

    cmd = f"""
        python3 generate_output_vcf.py --msp {msp} --tractor-output {tractor_output} {"--is-zipped" if zipped else ""} --output-path {v.ofile}
        """
    v.command(cmd)
    return v


def main(args):
    """Run batch LAI pipeline."""
    contig = args.contig
    logger.info(f"Running gnomAD LAI on chr{contig}")
    with run_batch(args, f"LAI - chr{contig}") as b:
        output_path = args.output_bucket

        if args.run_eagle:
            if args.sample_vcf:
                vcf = b.read_input(args.sample_vcf)
                e = eagle(b, vcf, contig, image=args.eagle_image)
                b.write_output(
                    e.ofile, dest=f"{output_path}eagle/output/phased_chr{contig}"
                )
            if args.ref_vcf:
                ref_vcf = b.read_input(args.ref_vcf)
                ref_e = eagle(b, ref_vcf, contig, image=args.eagle_image)
                b.write_output(
                    ref_e.ofile, dest=f"{output_path}eagle/phased_reference_chr{contig}"
                )

        if args.run_rfmix or args.run_xgmix:
            sample_map = b.read_input(args.pop_sample_map)
            genetic_map = b.read_input(args.genetic_map)
            phased_ref_vcf = (
                b.read_input(args.phased_ref_vcf)
                if args.phased_ref_vcf
                else ref_e.ofile.vcf
            )
            phased_sample_vcf = (
                b.read_input(args.phased_sample_vcf)
                if args.phased_sample_vcf
                else e.ofile.vcf
            )

            if args.run_rfmix:
                lai = rfmix(
                    b,
                    phased_sample_vcf,
                    phased_ref_vcf,
                    contig,
                    sample_map,
                    genetic_map,
                    image=args.rfmix_image,
                )
                b.write_output(lai.ofile, dest=f"{output_path}rfmix/output/chr{contig}")
            if args.run_xgmix:
                lai = xgmix(
                    b,
                    phased_sample_vcf,
                    genetic_map,
                    contig,
                    phased_ref_vcf,
                    sample_map,
                    threads=16,
                    mem="highmem",
                    image=args.xgmix_image,
                )
                b.write_output(lai.ofile, dest=f"{output_path}xgmix/output/chr{contig}")

        if args.run_tractor:
            phased_sample_vcf = (
                b.read_input_group(**{"vcf.gz": args.phased_sample_vcf})
                if args.phased_sample_vcf
                else e.ofile.vcf
            )
            msp_file = (
                b.read_input_group(**{"msp.tsv": args.msp_file})
                if args.msp_file
                else lai.ofile
            )
            t = tractor(
                b,
                msp_file,
                phased_sample_vcf,
                args.ancs,
                zipped=True,
                zip_output=args.zip_tractor_output,
                contig=contig,
                image=args.tractor_image,
            )
            b.write_output(t.ofile, dest=f"{output_path}tractor/output/chr{contig}")

        if args.make_lai_vcf:
            msp_file = b.read_input(args.msp_file) if args.msp_file else lai.ofile
            rg_def = {}
            if args.tractor_output:
                for i in range(args.ancs):
                    rg_def[
                        f"vcf{i}.gz"
                    ] = f"{args.tractor_output}.anc{i}.vcf{'.gz' if args.zip_tractor_output else ''}"
                    rg_def[
                        f"dos{i}.txt.gz"
                    ] = f"{args.tractor_output}.anc{i}.dosage.txt{'.gz' if args.zip_tractor_output else ''}"
                    rg_def[
                        f"ancdos{i}.txt.gz"
                    ] = f"{args.tractor_output}.anc{i}.hapcount.txt{'.gz' if args.zip_tractor_output else ''}"
            tractor_output = (
                b.read_input_group(**rg_def) if args.tractor_output else t.ofile
            )
            v = generate_lai_vcf(
                b,
                msp_file,
                tractor_output,
                zipped=args.zip_tractor_output,
                contig=contig,
            )
            b.write_output(
                v.ofile, dest=f"{output_path}tractor/output/chr{contig}_annotated"
            )


if __name__ == "__main__":
    p = init_arg_parser(
        default_cpu=8,
        default_billing_project="broad-mpg-gnomad",
        default_temp_bucket="gnomad-batch",
    )
    phasing_args = p.add_argument_group("Phasing", "Arguments for phasing samples")
    lai_args = p.add_argument_group(
        "Local Ancestry Inference",
        "Arguments for running local ancestry inference tools (rfmix, xgmix) on samples",
    )
    tractor_args = p.add_argument_group(
        "Tractor", "Arguments for running Tractor on samples"
    )
    vcf_args = p.add_argument_group("LAI VCF", "Arguments for generating LAI VCF")
    p.add_argument(
        "--output-bucket",
        required=True,
        help="Google bucket path for results. Each tool will create a subfolder here.",
    )
    phasing_args.add_argument(
        "--run-eagle",
        required=False,
        action="store_true",
        help="Whether to run eagle to phase samples.",
    )
    phasing_args.add_argument(
        "--sample-vcf",
        required=False,
        help="Google bucket path to sample VCF to phase.",
    )
    phasing_args.add_argument(
        "--ref-vcf",
        required=False,
        help="Google bucket path reference VCF to phase if separate.",
    )
    phasing_args.add_argument(
        "--eagle-image",
        help="Docker image for Eagle.",
        default="gcr.io/broad-mpg-gnomad/lai_phasing:latest",
    )
    lai_args.add_argument(
        "--phased-ref-vcf",
        required=False,
        help="Phased reference VCF, if supplied, will not re-run phasing.",
    )
    p.add_argument(
        "--phased-sample-vcf",
        help="VCF of phased samples, needed for LAI and/or Tractor runs.",
    )
    lai_args.add_argument(
        "--run-rfmix",
        required=False,
        action="store_true",
        help="Run local ancestry tool RFMix2.",
    )
    lai_args.add_argument(
        "--rfmix-image",
        help="Docker image for RFMix_v2.",
        default="gcr.io/broad-mpg-gnomad/lai_rfmix:latest",
    )
    lai_args.add_argument(
        "--run-xgmix",
        required=False,
        action="store_true",
        help="Run local ancestry tool XGMix.",
    )
    lai_args.add_argument(
        "--xgmix-image",
        help="Docker image for XGMix.",
        default="gcr.io/broad-mpg-gnomad/lai_xgmix:latest",
    )
    lai_args.add_argument(
        "--genetic-map",
        required=False,
        help="Genetic map required for RFMix_v2.",
        default="gs://gnomad-batch/mwilson/lai/inputs/rfmix/genetic_map_hg38.txt",
    )
    lai_args.add_argument(
        "--pop-sample-map", required=False, help="Sample population mapping for RFMix2."
    )
    p.add_argument("--contig", required=True, help="Chromosome to run LAI on.")
    p.add_argument(
        "--msp-file",
        required=False,
        help="Output from LAI program like RFMix_v2. Needed for Tractor and/or VCF generation.",
    )
    tractor_args.add_argument(
        "--run-tractor",
        required=False,
        action="store_true",
        help="Run Tractor's ExtractTracts.py script.",
    )
    tractor_args.add_argument(
        "--tractor-image",
        help="Docker image for Tractor.",
        default="gcr.io/broad-mpg-gnomad/lai_tractor:latest",
    )
    tractor_args.add_argument(
        "--ancs",
        required=False,
        help="Number of ancestries within the reference panel. Used to extract ancestry tracts from phased VCF in Tractor.",
        default=3,
        type=int,
    )
    tractor_args.add_argument(
        "--zip-tractor-output", help="Zip Tractors output", action="store_true"
    )
    vcf_args.add_argument(
        "--tractor-output",
        help="Path to tractor output files without anc.hapcount.txt and anc.dosage.txt, e.g. /Tractor/output/test_run",
    )
    vcf_args.add_argument(
        "--make-lai-vcf",
        help="Generate single VCF with ancestry AFs from tractor output, e.g. /Tractor/output/test_run",
        action="store_true",
    )
    p.add_argument(
        "--slack-channel",
        required=False,
        help="Slack channel to send job status to, needs @ for DM.",
    )
    args = p.parse_args()
    check_args(p, args)

    if args.slack_channel:
        from slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)

# noqa: D100
import argparse
import logging
from typing import Any

import hailtop.batch as hb

from tgg.batch.batch_utils import init_arg_parser, run_batch

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


def check_args(parser: argparse.ArgumentParser(), args: Any) -> None:
    """
    Check passed args to ensure pipeline can run properly.

    :param parser: Arg parser.
    :param args: Args from argparser.
    :return: None; will print error to stdout if arguments do not pass checks.
    """
    if not (
        args.run_eagle
        or args.split_phased_vcf
        or args.run_rfmix
        or args.run_xgmix
        or args.run_tractor
        or args.make_lai_vcf
    ):
        parser.error(
            "Need to specify at least one step to run (--run-eagle, --run-rfmix, --run-xgmix, --run-tractor, and/or --make-lai-vcf)."
        )
    if args.run_eagle:
        if not (args.cohort_vcf or args.ref_vcf):
            parser.error(
                "Need to specify either sample and/or reference vcfs (--cohort-vcf or --ref-vcf)."
            )
    if args.run_rfmix and args.run_xgmix:
        parser.error(
            "Can only specify one LAI tool, either RFMix or XGMix (--run-rfmix or --run-xgmix)."
        )
    if args.run_rfmix or args.run_xgmix:
        if not ((args.run_eagle and args.cohort_vcf) or args.phased_cohort_vcf):
            parser.error(
                "Need to specify either cohort vcf for eagle to run or pass a phased cohort vcf for RFMix to run (--run-eagle and --cohort-vcf or --phased-cohort-vcf)."
            )
        if not ((args.run_eagle and args.ref_vcf) or args.phased_ref_vcf):
            parser.error(
                "Need to specify either reference vcf for eagle to run or pass a phased cohort vcf for RFMix to run (--run-eagle and --reference-vcf or --phased-reference-vcf)."
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
        if not ((args.run_eagle and args.cohort_vcf) or args.phased_cohort_vcf):
            parser.error(
                "Need to specify either cohort vcf for eagle to run or pass a phased cohort vcf for RFMix to run (--run-eagle and --cohort-vcf or --phased-cohort-vcf)."
            )
        if not (args.run_rfmix or args.run_xgmix) and not args.msp_file:
            parser.error(
                "Need to run RFMix to generate MSP file or pass the MSP tsv file to the script, --msp-file."
            )
        if not args.n_ancs:
            parser.error(
                "Need to specify either number of continental ancestries within RFMix and phased cohort VCF."
            )


def split_vcf(
    batch: hb.Batch,
    phased_vcf: str,
    meta_table: str,
    contig: str,
    sample_type: str,
    mem: str = "highmem",
    storage: str = "100G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/tgg-methods-vm:latest",
) -> hb.job.Job:
    """
    Subset a VCF based on a provided sample list and process it using Hail Batch.

    :param meta_table: Path to the meta table for samples to subset to.
    :param phased_vcf: Path to the input VCF file.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param sample_type: Type of samples to subset to, i.e. cohort or reference.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "100G".
    :param cpu: The number of CPUs requested which is also used for threading, defaults to 16.
    :param image: Docker image for eagle job, defaults to "gcr.io/broad-mpg-gnomad/lai_phasing:latest".
    :return: Batch job.
    """
    split = batch.new_job(name="Run split_vcf")
    split.image(image)  # Set Docker image - I just put this as a placeholder
    split.memory(mem)  # Set memory requirement
    split.storage(storage)  # Set storage requirement
    split.cpu(cpu)  # Set CPU requirement
    split.declare_resource_group(
        ofile={"vcf.bgz": "{root}.vcf.bgz"}
    )  # using recode on this works if we arent piping to bgzip

    # Pipe to bgzip and index the VCF by vcftools
    # Command to execute the split_vcf function
    cmd = f"""vcftools --gzvcf {phased_vcf} \
        --keep {meta_table} \
        --recode \
        --stdout | bgzip -c > {split.ofile['vcf.bgz']}

        tabix -p vcf {split.ofile['vcf.bgz']}
        """
    split.command(cmd)

    return split


def eagle(
    batch: hb.Batch,
    vcf: str,
    contig: str,
    mem: str = "highmem",
    storage: str = "100G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/lai_phasing:latest",
) -> hb.Batch.new_job:
    """
    Run the phasing tool Eagle on passed VCF.

    :param batch: Hail batch object.
    :param vcf: VCF to phase.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "100G".
    :param cpu: The number of CPUs requested which is also used for threading, defaults to 16.
    :param image: Docker image for eagle job, defaults to "gcr.io/broad-mpg-gnomad/lai_phasing:latest".
    :return: Batch job.
    """
    e = batch.new_job(name=f"Eagle - chr{contig}")
    e.memory(mem)
    e.storage(storage)
    e.cpu(cpu)
    e.image(image)
    e.declare_resource_group(ofile={"vcf.gz": "{root}.vcf.gz"})

    cmd = f"""
    Eagle_v2.4.1/eagle \
        --geneticMapFile Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
        --numThreads {cpu} \
        --outPrefix {e.ofile} \
        --vcfOutFormat z \
        --vcf {vcf}
    """

    e.command(cmd)
    return e


def rfmix(
    batch: hb.Batch,
    cohort_pvcf: str,
    ref_pvcf: str,
    contig: str,
    sample_map: str,
    rf_genetic_map: str,
    mem: str = "highmem",
    storage: str = "100G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/lai_rfmix:latest",
) -> hb.Batch.new_job:
    """
    Run RFMix2 on phased VCF.

    :param batch: Hail batch object.
    :param cohort_pvcf: Phased cohort sample VCF from phasing tool like Eagle or SHAPEIT.
    :param ref_pvcf: Phased reference sample VCF from phasing tool like Eagle or SHAPEIT.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param sample_map: TSV file containing a mapping from sample IDs to ancestral populations, i.e. NA12878    EUR.
    :param rf_genetic_map: HapMap genetic map from SNP base pair positions to genetic coordinates in centimorgans.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "100G".
    :param cpu: The number of CPUs requested which is also used for threading, defaults to 16.
    :param image: RFMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_rfmix:latest".
    :return: Hail batch job.
    """
    r = batch.new_job(name=f"RFMix - chr{contig}")
    r.memory(mem)
    r.storage(storage)
    r.cpu(cpu)
    r.image(image)
    r.declare_resource_group(
        ofile={"msp.tsv": "{root}.msp.tsv", "fb.tsv": "{root}.fb.tsv"}
    )

    cmd = f"""
    ./rfmix \
        -f {cohort_pvcf} \
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
    cohort_pvcf: str,
    xg_genetic_map: str,
    contig: str,
    ref_pvcf: str,
    sample_map: str,
    mem: str = "highmem",
    storage: str = "100G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/lai_xgmix:latest",
) -> hb.Batch.new_job:
    """
    Run XGMix on phased VCF.

    :param batch: Hail batch object.
    :param cohort_pvcf: Phased cohort sample VCF from phasing tool like Eagle or SHAPEIT.
    :param xg_genetic_map: HapMap genetic map from SNP base pair positions to genetic coordinates in cM.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param ref_pvcf: Phased reference sample VCF from phasing tool like Eagle or SHAPEIT.
    :param sample_map: TSV file containing a mapping from sample IDs to ancestral populations, i.e. NA12878    EUR.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "100G".
    :param cpu: Number of CPUs requested, defaults to 16.
    :param image: XGMix Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_xgmix:latest".
    :return: Hail batch job.
    """
    x = batch.new_job(name=f"XGMix - chr{contig}")
    x.memory(mem)
    x.storage(storage)
    x.cpu(cpu)
    x.image(image)
    x.declare_resource_group(
        ofile={"msp.tsv": "{root}.msp.tsv", "fb.tsv": "{root}.fb.tsv"}
    )

    cmd = f"""
    python3 XGMIX.py {cohort_pvcf} {xg_genetic_map} /io/tmp/xgmix/output chr{contig} False {ref_pvcf} {sample_map}
    ln -s /io/tmp/xgmix/output/output.msp.tsv {x.ofile['msp.tsv']}
    ln -s /io/tmp/xgmix/output/output.fb.tsv {x.ofile['fb.tsv']}
    """

    x.command(cmd)
    return x


def tractor(
    batch: hb.Batch,
    msp: str,
    cohort_pvcf: str,
    n_ancs: int,
    input_zipped: bool,
    zip_output: bool,
    contig: str,
    mem: str = "highmem",
    storage: str = "200G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/lai_tractor:latest",
) -> hb.Batch.new_job:
    """
    Run Tractor's ExtractTract.py script.

    :param batch: Hail batch object.
    :param msp: MSP tsv file from LAI tool like RFMix2 or XGMix.
    :param vcf: Phased cohort sample VCF from phasing tool like Eagle or SHAPEIT.
    :param n_ancs: Number of ancestral populations within the MSP file.
    :param input_zipped: Whether the input VCF file is zipped or not, i.e. ends in vcf.gz.
    :param zip_output: Whether to zip the tool's output files.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "200G".
    :param cpu: The number of CPUs requested which is also used for threading, defaults to 16.
    :param image: Tractor Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_tractor:latest".
    :return: Hail Batch job.
    """
    t = batch.new_job(name=f"Tractor - chr{contig}")
    t.memory(mem)
    t.storage(storage)
    t.cpu(cpu)
    t.image(image)
    rg_def = {}
    file_extension = ".gz" if zip_output else ""
    for i in range(n_ancs):
        rg_def[f"vcf{i}{file_extension}"] = f"{{root}}.anc{i}.vcf{file_extension}"
        rg_def[
            f"dos{i}.txt{file_extension}"
        ] = f"{{root}}.anc{i}.dosage.txt{file_extension}"
        rg_def[
            f"ancdos{i}.txt{file_extension}"
        ] = f"{{root}}.anc{i}.hapcount.txt{file_extension}"

    t.declare_resource_group(ofile=rg_def)
    input_zipped = "--zipped" if input_zipped else ""
    zip_output = "--zip-output" if zip_output else ""

    cmd = f"""
    python3 ExtractTracts.py --msp {msp} --vcf {cohort_pvcf} --num-ancs={n_ancs} {input_zipped} {zip_output} --output-path={t.ofile}
    """

    t.command(cmd)
    return t


def generate_lai_vcf(
    batch: hb.Batch,
    msp: str,
    tractor_output: str,
    input_zipped: bool,
    contig: str,
    mt_path_for_adj: str,
    add_gnomad_af: bool,
    mem: str = "highmem",
    storage: str = "200G",
    cpu: int = 16,
    image: str = "gcr.io/broad-mpg-gnomad/lai_vcf:latest",
) -> hb.Batch.new_job:
    """
    Run generate_output_vcf.py script.

    :param batch: Hail batch object.
    :param msp: MSP tsv file from LAI tool like RFMix2 or XGMix.
    :param tractor_output: Path to Tractor's output files.
    :param input_zipped: Whether the input VCF file is zipped or not, i.e. ends in vcf.gz.
    :param contig: Which chromosome the VCF contains. This must be a single chromosome.
    :param mt_path_for_adj: Path to MT to filter to high quality genotypes before calculating AC.
    :param add_gnomad_af: Whether to add gnomAD's population AFs for AMR, NFE, AFR, and EAS.
    :param mem: Hail batch job memory, defaults to "highmem".
    :param storage: Hail batch job storage, defaults to "200G".
    :param cpu: The number of CPUs requested which is also used for threading, defaults to 16.
    :param image: VCF Docker image, defaults to "gcr.io/broad-mpg-gnomad/lai_vcf:latest".
    :return: Hail Batch job.
    """
    v = batch.new_job(name=f"Generate final VCF - chr{contig}")
    v.memory(mem)
    v.storage(storage)
    v.cpu(cpu)
    v.image(image)
    v.declare_resource_group(ofile={"vcf.bgz": "{root}_lai_annotated.vcf.bgz"})

    if mt_path_for_adj:
        mt_path_for_adj = f"--mt-path-for-adj {mt_path_for_adj}"

    cmd = f"""
    python3 generate_output_vcf.py --msp {msp} --tractor-output {tractor_output} {"--is-zipped" if input_zipped else ""} {mt_path_for_adj if mt_path_for_adj else ""} {"--add-gnomad-af" if add_gnomad_af else ""} --output-path {v.ofile}
    """

    v.command(cmd)
    return v


def main(args):
    """Run batch local ancestry inference (LAI) pipeline.

    The pipeline has four steps that can run independently or in series using either user input or a previous step's output:

        - Phase a cohort VCF and a reference VCF using Eagle.
        - Run a local ancestry tool, either RFMix or XGMix, on phased cohort VCF.
        - Run Tractor to extract ancestral components from the phased cohort VCF and generate a VCF, dosage counts, and haplotype counts per ancestry.
        - Generate a single VCF with ancestry-specific call statistics (AC, AN, AF).
    """
    contig = args.contig
    contig = contig[3:] if contig.startswith("chr") else contig
    logger.info("Running gnomAD LAI on chr%s", contig)
    with run_batch(args, f"LAI - chr{contig}") as b:
        output_path = args.output_bucket

        if args.run_eagle:
            if args.cohort_vcf:
                logger.info("Running eagle on cohort VCF...")
                vcf = b.read_input(args.cohort_vcf)
                e = eagle(
                    b,
                    vcf,
                    contig,
                    mem=args.eagle_mem,
                    storage=args.eagle_storage,
                    cpu=args.eagle_cpu,
                    image=args.eagle_image,
                )
                b.write_output(
                    e.ofile,
                    dest=f"{output_path}chr{contig}/eagle/output/phased_chr{contig}",
                )
            if args.ref_vcf:
                logger.info("Running eagle on reference VCF...")
                ref_vcf = b.read_input(args.ref_vcf)
                ref_e = eagle(
                    b,
                    ref_vcf,
                    contig,
                    mem=args.eagle_mem,
                    storage=args.eagle_storage,
                    cpu=args.eagle_cpu,
                    image=args.eagle_image,
                )
                b.write_output(
                    ref_e.ofile,
                    dest=f"{output_path}chr{contig}/eagle/phased_reference_chr{contig}",
                )

        if args.split_phased_vcf:
            logger.info("Splitting phased VCF...")
            # Define lists to store sample IDs for cohort and reference panels
            ref_meta_table = b.read_input(args.ref_meta_for_split)
            cohort_meta_table = b.read_input(args.cohort_meta_for_split)

            phased_vcf = (
                b.read_input(args.phased_cohort_vcf)
                if args.phased_cohort_vcf
                else e.ofile["vcf.gz"]
            )

            split_cohort_vcf = split_vcf(
                b, phased_vcf, cohort_meta_table, contig=contig, sample_type="cohort"
            )
            b.write_output(
                split_cohort_vcf.ofile,
                dest=f"{output_path}chr{contig}/split_phased_vcf/cohort",
            )

            split_ref_vcf = split_vcf(
                b, phased_vcf, ref_meta_table, contig=contig, sample_type="reference"
            )
            b.write_output(
                split_ref_vcf.ofile,
                dest=f"{output_path}chr{contig}/split_phased_vcf/reference",
            )

        if args.run_rfmix or args.run_xgmix:
            sample_map = b.read_input(args.pop_sample_map)
            genetic_map = b.read_input(args.genetic_map)
            phased_ref_vcf = (
                b.read_input(args.phased_ref_vcf)
                if args.phased_ref_vcf
                else ref_e.ofile["vcf.gz"]
            )
            phased_cohort_vcf = (
                b.read_input(args.phased_cohort_vcf)
                if args.phased_cohort_vcf
                else e.ofile["vcf.gz"]
            )

            if args.run_rfmix:
                logger.info("Running Local Ancestry Inference tool RFMix v2...")
                lai = rfmix(
                    b,
                    phased_cohort_vcf,
                    phased_ref_vcf,
                    contig,
                    sample_map,
                    genetic_map,
                    mem=args.lai_mem,
                    storage=args.lai_storage,
                    cpu=args.lai_cpu,
                    image=args.rfmix_image,
                )
                b.write_output(
                    lai.ofile, dest=f"{output_path}chr{contig}/rfmix/output/chr{contig}"
                )
            if args.run_xgmix:
                logger.info("Running Local Ancestry Inference tool XGMix...")
                lai = xgmix(
                    b,
                    phased_cohort_vcf,
                    genetic_map,
                    contig,
                    phased_ref_vcf,
                    sample_map,
                    mem=args.lai_mem,
                    storage=args.lai_storage,
                    cpu=args.lai_cpu,
                    image=args.xgmix_image,
                )
                b.write_output(
                    lai.ofile, dest=f"{output_path}chr{contig}/xgmix/output/chr{contig}"
                )

        if args.run_tractor:
            logger.info("Running Tractor...")
            # Both inputs have a specified extension so batch can find the file and pass it to Tractor which expects files without extensions
            msp_file = (
                b.read_input_group(**{"msp.tsv": args.msp_file})
                if args.msp_file
                else lai.ofile
            )
            phased_cohort_vcf = (
                b.read_input_group(**{"vcf.gz": args.phased_cohort_vcf})
                if args.phased_cohort_vcf
                else e.ofile
            )
            t = tractor(
                b,
                msp_file,
                phased_cohort_vcf,
                args.n_ancs,
                input_zipped=True,
                zip_output=args.zip_tractor_output,
                contig=contig,
                mem=args.tractor_mem,
                storage=args.tractor_storage,
                cpu=args.tractor_cpu,
                image=args.tractor_image,
            )
            b.write_output(
                t.ofile, dest=f"{output_path}chr{contig}/tractor/output/chr{contig}"
            )

        if args.make_lai_vcf:
            logger.info("Generating output VCF...")
            msp_file = (
                b.read_input(args.msp_file) if args.msp_file else lai.ofile["msp.tsv"]
            )
            rg_def = {}
            if args.tractor_output:
                for i in range(args.n_ancs):
                    rg_def[
                        f"anc{i}.dosage.txt"
                    ] = f"{args.tractor_output}.dos{i}.txt{'.gz' if args.zip_tractor_output else ''}"
                    rg_def[
                        f"anc{i}.hapcount.txt"
                    ] = f"{args.tractor_output}.ancdos{i}.txt{'.gz' if args.zip_tractor_output else ''}"
            tractor_output = (
                b.read_input_group(**rg_def) if args.tractor_output else t.ofile
            )
            v = generate_lai_vcf(
                b,
                msp_file,
                tractor_output,
                input_zipped=args.zip_tractor_output,
                contig=contig,
                mt_path_for_adj=args.mt_path_for_adj,
                add_gnomad_af=args.add_gnomad_af,
                mem=args.vcf_mem,
                storage=args.vcf_storage,
                cpu=args.vcf_cpu,
                image=args.vcf_image,
            )
            b.write_output(
                v.ofile,
                dest=f"{output_path}chr{contig}/tractor/output/chr{contig}_annotated{'_adj'if args.mt_path_for_adj else ''}",
            )
    logger.info("Batch LAI pipeline run complete!")


if __name__ == "__main__":
    p = init_arg_parser(
        default_cpu=16,
        default_billing_project="gnomad-production",
        default_temp_bucket="gnomad-batch",
    )
    multi_args = p.add_argument_group(
        "Multi-step use", "Arguments used by multiple steps"
    )
    multi_args.add_argument(
        "--contig",
        required=True,
        help="Chromosome to run LAI on with the 'chr' prefix.",
    )
    multi_args.add_argument(
        "--output-bucket",
        required=True,
        help="Google bucket path with final / included. Each steps' result will be written to within a chromosome subfolder here.",
    )
    multi_args.add_argument(
        "--phased-cohort-vcf",
        required=False,
        help="Zipped VCF of phased cohort samples, needed for LAI and/or Tractor runs.",
    )
    multi_args.add_argument(
        "--msp-file",
        required=False,
        help="Output from LAI program like RFMix_v2. Needed for Tractor and/or VCF generation.",
    )
    phasing_args = p.add_argument_group("Phasing", "Arguments for phasing samples")
    phasing_args.add_argument(
        "--run-eagle",
        required=False,
        action="store_true",
        help="Whether to run eagle to phase samples.",
    )
    phasing_args.add_argument(
        "--eagle-mem",
        default="highmem",
        help="Memory for eagle batch job.",
    )
    phasing_args.add_argument(
        "--eagle-storage",
        default="100G",
        help="Storage for eagle batch job.",
    )
    phasing_args.add_argument(
        "--eagle-cpu",
        default=16,
        help="CPU for eagle batch job.",
    )
    phasing_args.add_argument(
        "--cohort-vcf",
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
    splitting_args = p.add_argument_group(
        "Splitting", "Arguments for splitting phased VCF by cohort and reference"
    )
    splitting_args.add_argument(
        "--ref-meta-for-split",
        required=False,
        help="TSV of reference samples to split phased VCF.",
    )
    splitting_args.add_argument(
        "--cohort-meta-for-split",
        required=False,
        help="TSV of cohort samples to split phased VCF.",
    )
    splitting_args.add_argument(
        "--split-phased-vcf",
        required=False,
        action="store_true",
        help="Whether to split phased VCF by cohort and reference.",
    )
    lai_args = p.add_argument_group(
        "Local Ancestry Inference",
        "Arguments for running local ancestry inference tools (rfmix, xgmix) on samples",
    )
    lai_args.add_argument(
        "--lai-mem",
        default="highmem",
        help="Memory for LAI tool batch job.",
    )
    lai_args.add_argument(
        "--lai-storage",
        default="100G",
        help="Storage for LAI tool batch job.",
    )
    lai_args.add_argument(
        "--lai-cpu",
        default=16,
        help="CPU for LAI tool batch job.",
    )
    lai_args.add_argument(
        "--phased-ref-vcf",
        required=False,
        help="Zipped VCF of phased reference samples. If supplied, the phasing step will not run on reference data.",
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
        help="Genetic map from SNP base pair positions to genetic coordinates in centimorgans, required for RFMix_v2.",
        default="gs://gnomad-batch/mwilson/lai/inputs/rfmix/genetic_map_hg38.txt",
    )
    lai_args.add_argument(
        "--pop-sample-map", required=False, help="Sample population mapping for RFMix2."
    )
    tractor_args = p.add_argument_group(
        "Tractor", "Arguments for running Tractor on samples"
    )
    tractor_args.add_argument(
        "--run-tractor",
        required=False,
        action="store_true",
        help="Run Tractor's ExtractTracts.py script.",
    )
    tractor_args.add_argument(
        "--tractor-mem",
        default="highmem",
        help="Memory for Tractor batch job.",
    )
    tractor_args.add_argument(
        "--tractor-storage",
        default="200G",
        help="Storage for Tractor batch job.",
    )
    tractor_args.add_argument(
        "--tractor-cpu",
        default=16,
        help="CPU for Tractor batch job.",
    )
    tractor_args.add_argument(
        "--tractor-image",
        help="Docker image for Tractor.",
        default="gcr.io/broad-mpg-gnomad/lai_tractor:latest",
    )
    tractor_args.add_argument(
        "--n-ancs",
        help="Number of ancestries within the reference panel. Used to extract ancestry tracts from phased VCF in Tractor.",
        default=3,
        type=int,
    )
    tractor_args.add_argument(
        "--zip-tractor-output", help="Zip Tractors output", action="store_true"
    )
    vcf_args = p.add_argument_group("LAI VCF", "Arguments for generating LAI VCF")
    vcf_args.add_argument(
        "--tractor-output",
        help="Path to tractor output files without anc.hapcount.txt and anc.dosage.txt, e.g. /Tractor/output/test_run",
    )
    vcf_args.add_argument(
        "--make-lai-vcf",
        help="Generate single VCF with ancestry AFs from tractor output.",
        action="store_true",
    )
    vcf_args.add_argument(
        "--mt-path-for-adj",
        help="Path to hail MatrixTable generated from pipeline input VCF. Must contain GT, GQ, DP, and AB fields. If MT path provided, script will filter to high quality GTs only.",
    )
    vcf_args.add_argument(
        "--add-gnomad-af",
        help="Add gnomAD population allele frequencies from AMR, AFR, EAS, and NFE  to output VCF.",
        action="store_true",
    )
    vcf_args.add_argument(
        "--vcf-mem",
        default="highmem",
        help="Memory for VCF generation batch job.",
    )
    vcf_args.add_argument(
        "--vcf-storage",
        default="200G",
        help="Storage for VCF generation batch job.",
    )
    vcf_args.add_argument(
        "--vcf-cpu",
        default=16,
        help="CPU for VCF generation batch job.",
    )
    vcf_args.add_argument(
        "--vcf-image",
        help="Docker image for VCF generation.",
        default="gcr.io/broad-mpg-gnomad/lai_vcf:latest",
    )
    args = p.parse_args()
    check_args(p, args)

    main(args)

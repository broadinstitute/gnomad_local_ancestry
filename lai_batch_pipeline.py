import logging
import hailtop.batch as hb

from gnomad.utils.slack import slack_notifications
from slack_creds import slack_token
from batch.batch_utils import (
    init_arg_parser,
    init_job,
    print_memory_stats,
    run_batch,
)

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

GCLOUD_USER_ACCOUNT = "mwilson@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://seqr-mwilson-misc/"
DOCKER_IMAGE = "gcr.io/broad-mpg-gnomad/lai_phasing:latest"


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
    :param contig: Which chromosome the VCF contains
    :param threads: The number of threads, should match the number of CPUs requested, defaults to 16
    :return: Batch job
    """
    e = batch.new_job(name=f"Eagle-chr{contig}")
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


def run_lai(args):
    contig = args.contig
    logger.info(f"Running gnomAD LAI on chr{contig}")
    with run_batch(args, f"LAI - chr{contig}") as b:
        vcf = b.read_input(args.sample_vcf)
        output_path = args.output_bucket

        if args.ref_vcf:
            ref_vcf = b.read_input(args.ref_vcf)
            ref_e = eagle(b, ref_vcf, contig)
            b.write_output(ref_e.ofile, dest=f"{output_path}eagle/chr{contig}_reference")

        e = eagle(b, vcf, contig)
        b.write_output(e.ofile, dest=f"{output_path}eagle/chr{contig}_amr")


if __name__ == "__main__":
    p = init_arg_parser(
        default_cpu=8,
        default_billing_project="broad-mpg-gnomad",
        default_temp_bucket="gnomad-batch",
    )
    p.add_argument(
        "--sample-vcf", required=True, help="Google bucket path to sample VCF to phase"
    )
    p.add_argument(
        "--ref-vcf",
        required=False,
        help="Google bucket path reference VCF to phase if separate",
    )
    p.add_argument(
        "--output-bucket", required=True, help="Google bucket path for results"
    )
    p.add_argument("--rf-genetic-map", required=False, help="Genetic map for required for RFMix2")
    p.add_argument("--contig", required=True, help="Chromosomes to run LAI on")
    p.add_argument("--slack-channel", help="Slack channel to send job status to, needs @ for DM")
    args = p.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            run_lai(args)
    else:
        run_lai(args)

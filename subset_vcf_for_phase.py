# noqa: D100
import logging

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.grch38.gnomad import public_release
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.pipeline import filter_rows_for_qc
from gnomad.utils.filtering import subset_samples_and_variants
from gnomad.utils.reference_genome import get_reference_genome
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)


def main(
    mt_path: str,
    samples_path: str,
    output_bucket: str,
    contigs: list,
    dense: bool = False,
    gt_expr: str = "LGT",
    min_callrate: float = 0.9,
    min_af: float = 0.001,
) -> hl.MatrixTable:
    """
    Subset a matrix table to specified samples and across specified contigs.

    Subset and filter on min_callrate of 0.9 and min_af of 0.001. Export each subsetted contig individually.

    :param mt_path: Path to MatrixTable to subset from.
    :param samples_path: Path to TSV of sample IDs to subset to. The TSV must have a header of 's'.
    :param output_bucket: Path to output bucket for contig MT and VCF.
    :param contigs: List of contigs as integers.
    :param dense: Boolean of whether source MT is dense. Defaults to False.
    :param gt_expr: Boolean of GT expression in MT. Defaults to 'LGT'.
    :param min_callrate: Minimum variant callrate for variant QC. Defaults to 0.9.
    :param min_af: Minimum allele frequency for variant QC. Defaults to 0.001.
    """
    logger.info("Running script on %s...", contigs)
    full_mt = hl.read_matrix_table(mt_path)
    for contig in contigs:
        contig = f"chr{contig}"
        logger.info("Subsetting %s...", contig)
        mt = hl.filter_intervals(
            full_mt,
            [
                hl.parse_locus_interval(
                    contig, reference_genome=get_reference_genome(full_mt.locus)
                )
            ],
        )
        if "s" not in hl.import_table(samples_path, no_header=False).row.keys():
            raise DataException(
                "The TSV provided by `sample_path` must include a header with a column labeled `s` for the sample IDs to keep in the subset."
            )

        mt = subset_samples_and_variants(
            mt, sample_path=samples_path, sparse=not dense, gt_expr=gt_expr
        )

        if not dense:
            mt = hl.MatrixTable(
                hl.ir.MatrixKeyRowsBy(
                    mt._mir, ["locus", "alleles"], is_sorted=True
                )  # Prevents hail from running sort on genotype MT which is already sorted by a unique locus
            )
            mt = mt.drop("gvcf_info")
            mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
            mt = hl.experimental.densify(mt)
            mt = mt.filter_rows(
                hl.len(mt.alleles) > 1
            )  # Note: This step is sparse-specific, removing monoallelic sites after densifying
        else:
            mt = hl.split_multi_hts(mt)

        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        logger.info(
            "Filtering to variants with greater than %d callrate and %d allele frequency",
            min_callrate,
            min_af,
        )
        if args.gnomad_release_only:
            logger.info("Filtering to gnomAD v3.1.1 release variants")
            mt = mt.filter_rows(
                hl.is_defined(public_release("genomes").ht()[mt.row_key])
            )
        mt = filter_rows_for_qc(
            mt,
            min_callrate=min_callrate,
            min_af=min_af,
            min_inbreeding_coeff_threshold=None,
            min_hardy_weinberg_threshold=None,
        )
        mt = mt.checkpoint(
            f"{output_bucket}{contig}/{contig}_dense_bia_snps.mt", overwrite=True,
        )
        logger.info(
            "Subsetted %s to %d variants and %d samples",
            contig,
            mt.rows().count(),
            mt.cols().count(),
        )
        hl.export_vcf(mt, f"{output_bucket}{contig}/{contig}_dense_bia_snps.vcf.bgz")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt-path", help="MatrixTable to subset from", required=True)
    parser.add_argument(
        "--samples-path",
        help="TSV of samples, expects the TSV to have a header with the label `s`",
    )
    parser.add_argument(
        "--output-bucket", help="Bucket for MTs and VCFs", required=True
    )
    parser.add_argument(
        "--dense", help="Whether MT is dense. Defaults to False", action="store_true"
    )
    parser.add_argument(
        "--gt-expr",
        help="Genotype expression, typically 'LGT' is for sparse MTs while 'GT' for dense.",
        default="LGT",
    )
    parser.add_argument(
        "--min-callrate",
        help="Minimum callrate threshiold as float for variant QC",
        type=float,
        default=0.90,
    )
    parser.add_argument(
        "--min-af",
        help="Minimum allele frequency as float for variant QC",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        help="Integer contigs to run subsetting on",
        required=True,
    )
    parser.add_argument(
        "--gnomad-release-only",
        help="Filter to only variants in the gnomad v3.1 release",
        action="store_true",
    )
    args = parser.parse_args()
    main(
        mt_path=args.mt_path,
        samples_path=args.samples_path,
        output_bucket=args.output_bucket,
        contigs=args.contigs,
        dense=args.dense,
        gt_expr=args.gt_expr,
        min_callrate=args.min_callrate,
        min_af=args.min_af,
    )

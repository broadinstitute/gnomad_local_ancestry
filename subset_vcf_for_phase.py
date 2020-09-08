import logging

import hail as hl

from gnomad.sample_qc.pipeline import filter_rows_for_qc
from gnomad.utils.filtering import subset_samples_and_variants

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def subset_and_export_chr_vcfs(
    mt_path: str,
    samples_path: str,
    output_bucket: str,
    contigs: list,
    sparse: bool = True,
    gt_expr: str = "LGT",
) -> hl.MatrixTable:
    """
    Subset a matrix table to specified samples and across specified contigs, filter
    on min_callrate of 0.9 and min_af of 0.001, then export each subsetted contig individually
    :param mt: Path to atrixTable to subset from.
    :param samples_path: Path to tsv of sample IDs with header "s".
    :param output_bucket: Path to output bucket for contig MT and VCF.
    :param contigs: List of contigs as integers.
    :param sparse: Boolean of whether source MT is sparse. Defaults to True.
    :param gt_expr: Boolean of GT expression in MT. Defaults to "LGT".
    """
    logger.info(f"Running script on {contigs}...")
    mt = hl.read_matrix_table(mt_path)
    for contig in contigs:
        contig = f"chr{contig}"
        logger.info(f"Subsetting {contig}...")
        mt = hl.filter_intervals(
            mt, [hl.parse_locus_interval(contig, reference_genome="GRCh38")]
        )
        mt = subset_samples_and_variants(
            mt, sample_path=samples_path, sparse=sparse, gt_expr=gt_expr
        )

        if sparse:
            mt = mt.key_rows_by('locus', 'alleles')
            mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
            mt = hl.experimental.densify(mt)
            mt = mt.filter_rows(hl.len(mt.alleles) > 1).drop(
                "gvcf_info"
            )  # Note: This step is sparse-specific, removing monoallelic sites after densifying
        else:
            mt = hl.split_multi_hts(mt)

        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        mt = filter_rows_for_qc(mt, min_callrate=0.9, min_af=0.001)

        mt = mt.checkpoint(
            f"{output_bucket}{contig}/gnomad_{contig}_dense_bia_snps.mt",  _read_if_exists=True
        )
        logger.info(f"Subsetted {contig} to {mt.rows().count()} variants and {mt.cols().count()} samples")
        hl.export_vcf(
            mt, f"{output_bucket}{contig}/gnomad_{contig}_dense_bia_snps.vcf.bgz"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt-path", help="MatrixTable to subset from", required=True)
    parser.add_argument("--samples-path", help="TSV of samples")
    parser.add_argument("--output-bucket", help="Bucket for MTs and VCFs", required=True)
    parser.add_argument("--sparse", help="Whether MT is sparse. Defaults to True", action="store_true")
    parser.add_argument("--gt-expr", help="Genotype expression, typically 'LGT' is for sparse MTs while 'GT' for dense.", default="LGT")
    parser.add_argument(
        "--contigs", nargs="+", help="Integer contigs to run subsetting on", required=True
    )
    args = parser.parse_args()
    subset_and_export_chr_vcfs(
        mt_path=args.mt_path,
        samples_path=args.samples_path,
        output_bucket=args.output_bucket,
        contigs=args.contigs,
        sparse=args.sparse,
        gt_expr=args.gt_expr,
    )

import logging as logger

import hail as hl

from gnomad.sample_qc.pipeline import filter_rows_for_qc
from gnomad.utils.filtering import subset_samples_and_variants


def subset_and_export_chr_vcfs(
    mt: hl.MatrixTable,
    samples_path: str,
    output_bucket: str,
    contigs: list,
    sparse: bool = True,
    gt_expr: str = "LGT",
) -> hl.MatrixTable:
    """
    Subset a matrix table to specified samples and across specified contigs, filter
    on min_callrate of 0.9 and min_af of 0.001, then export each subsetted contig individually
    :param mt: MatrixTable to subset from.
    :param samples_path: Path to tsv of sample IDs with header "s".
    :param output_bucket: Path to output bucket for contig MT and VCF.
    :param contigs: List of contigs as integers.
    :param sparse: Boolean of whether source MT is sparse. Defaults to True.
    :param gt_expr: Boolean of GT expression in MT. Defaults to "LGT".
    """
    for contig in contigs:
        contig = f"chr{contig}"
        mt = hl.filter_intervals(
            mt, [hl.parse_locus_interval(contig, reference_genome="GRCh38")]
        )
        mt = subset_samples_and_variants(
            mt, samples_path, sparse=sparse, gt_expr=gt_expr
        )

        if sparse:
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
            f"{output_bucket}{contig}/gnomad_{contig}_dense_bia_snps.mt", overwrite=True
        )
        hl.export_vcf(
            mt, f"{output_bucket}{contig}/gnomad_{contig}_dense_bia_snps.vcf.bgz"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt", help="MatrixTable to subset from")
    parser.add_argument("--samples-path", help="TSV of samples")
    parser.add_argument("--output-bucket", help="Bucket for MTs and VCFs")
    parser.add_argument("--sparse", help="Whether MT is sparse")
    parser.add_argument("--gt-expr", help="Ancestral pops to subset")
    parser.add_argument(
        "--contigs", nargs="?", help="comma separated integer contigs to run"
    )
    args = parser.parse_args()
    subset_and_export_chr_vcfs(
        args.mt,
        args.samples_path,
        args.output_bucket,
        args.contigs,
        args.sparse,
        args.gt_expr,
    )

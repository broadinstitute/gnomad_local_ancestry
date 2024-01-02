# noqa: D100
import logging
from typing import List, Optional

import hail as hl
from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.sample_qc.pipeline import filter_rows_for_qc


from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
from gnomad_qc.v4.resources.meta import meta

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)

CONTIGS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
]


def get_subset_samples(
    samples_path: Optional[str] = None,
    pops: Optional[List[str]] = None,
    hgdp: bool = False,
    tgp: bool = False,
) -> hl.Table:
    """
    Get samples to subset out of v3 VDS.

    :param samples_path: Path to TSV of sample IDs to subset to. The TSV must have a header of 's'.
    :param pops: Optional list of populations to include in subset. Defaults to None.
    :param hgdp: Boolean of whether to include HGDP in subset. Defaults to False.
    :param tgp: Boolean of whether to include TGP in subset. Defaults to False.
    """
    meta_ht = meta(data_type="genomes").ht()
    samples_to_keep = []

    def _add_filtered_meta(condition):
        ht_to_add = meta_ht.filter(condition).select()
        samples_to_keep.append(ht_to_add)

    if samples_path:
        sample_ht = hl.import_table(samples_path)
        sample_ht = sample_ht.select()

    # Select released samples with infered population
    if pops:
        pops_to_keep = hl.literal(pops)
        _add_filtered_meta(
            (pops_to_keep.contains(meta_ht.population_inference.pop))
            & (meta_ht.release)
        )

    # Keep HGDP samples in subset
    if hgdp:
        _add_filtered_meta(meta_ht.subsets.hgdp)

    # Keep TGP samples in subset
    if tgp:
        _add_filtered_meta(meta_ht.subsets.tgp)

    # Create final subset sample table
    sample_ht = hl.Table.union(*samples_to_keep)
    return sample_ht


def main(args):
    """
    Subset a matrix table to specified samples and across specified contigs.

    Subset and filter on min_callrate of 0.9 and min_af of 0.001. Export each subsetted contig individually.

    :param mt_path: Path to MatrixTable to subset from.
    :param samples_path: Path to TSV of sample IDs to subset to. The TSV must have a header of 's'.
    :param output_bucket: Path to output bucket for contig MT and VCF.
    :param contigs: List of contigs as integers.
    :param dense: Boolean of whether source MT is dense. Defaults to False.
    :param min_callrate: Minimum variant callrate for variant QC. Defaults to 0.9.
    :param min_af: Minimum allele frequency for variant QC. Defaults to 0.001.
    """
    hl.init(
        log="/subset_vcf_for_phase.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day/subset_vcf_for_phase/",
    )

    output_path = args.output_path
    contigs = args.contigs
    min_callrate = args.min_callrate
    min_af = args.min_af
    test = args.test
    contigs = CONTIGS if not contigs else contigs

    # TODO: Come back to args needed in get_gnomad_v3_vds
    logger.info("Running script on %s...", contigs)
    vds = get_gnomad_v3_vds()

    if test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(2)),
            vds.variant_data._filter_partitions(range(2)),
        )

    sample_ht = get_subset_samples(
        samples_path=args.samples_path, pops=args.pops, hgdp=args.hgdp, tgp=args.tgp
    )

    logger.info("Subsetting to %d samples", sample_ht.count())
    vds = hl.vds.filter_samples(vds, sample_ht, remove_dead_alleles=True)

    logger.info(
        "Applying min_rep to the variant data MT because remove_dead_alleles may "
        "result in variants that do not have the minimum representation."
    )
    vd = vds.variant_data
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vd.key_rows_by(**hl.min_rep(vd.locus, vd.alleles)),
    )

    for contig in contigs:
        logger.info("Subsetting %s...", contig)
        vds = hl.vds.filter_intervals(
            vds,
            [hl.parse_locus_interval(contig, reference_genome="GRCh38")],
        )
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.drop("gvcf_info")

        # TODO: See if this is needed: mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        logger.info(
            "Filtering to variants with greater than %f callrate and %f allele frequency",
            min_callrate,
            min_af,
        )
        mt = filter_rows_for_qc(
            mt,
            min_callrate=min_callrate,
            min_af=min_af,
            min_inbreeding_coeff_threshold=None,
            min_hardy_weinberg_threshold=None,
        )
        mt = mt.checkpoint(
            f"{output_path}/{contig}/{contig}_dense_bia_snps.mt",
            overwrite=True,
        )

        logger.info("Exporting VCF for %s...", contig)
        hl.export_vcf(
            mt, f"{output_path}{contig}/{contig}_dense_bia_snps.vcf.bgz", tabix=True
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt-path", help="MatrixTable to subset from", required=True)
    parser.add_argument(
        "--samples-path",
        help="TSV of samples, expects the TSV to have a header with the label `s`",
    )
    parser.add_argument("--output-path", help="Bucket for MTs and VCFs", required=True)
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
    )
    parser.add_argument(
        "--pops",
        nargs="+",
        help="Populations to include in subset",
    )
    parser.add_argument(
        "--tgp",
        help="Include TGP in subset",
        action="store_true",
    )
    parser.add_argument(
        "--hgdp",
        help="Include HGDP in subset",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Subset to 2 partitions",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)

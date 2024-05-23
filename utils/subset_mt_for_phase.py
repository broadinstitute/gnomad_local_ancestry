"""Module to subset MT for phasing."""

import logging
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

MT_PATH = "gs://gnomad-lai/afr/sample_subsets/{contig}/{contig}_dense_biallelic_snps.mt"
META_PATH = "gs://gnomad-lai/afr/sample_subsets/afr_subset_samples_with_hgdp_tgp.ht"


def subset_mt(mt: hl.MatrixTable, populations: list, output_path: str):
    """
    Subset a MatrixTable based on specified populations and create a sample map.

    :param mt: Input MatrixTable.
    :param populations: List of populations to subset (e.g., ["afr", "nfe"]).
    :param output_file_path: Path for the output sample map file.
    :return: Filtered MatrixTable.
    """
    logger.info("Reading in meta for subsetting...")
    meta = hl.read_table(META_PATH)

    logger.info("Annotating MatrixTable with meta information...")
    mt = mt.annotate_cols(
        hgdp=meta[mt.s].subsets.hgdp,
        tgp=meta[mt.s].subsets.tgp,
        project_pop=meta[mt.s].project_meta.project_pop,
        project_subpop=meta[mt.s].project_meta.project_subpop,
    )
    logger.info(f"{mt.count_cols()} samples in MatrixTable.")
    logger.info("Removing HGDP/TGP samples from unspecified populations...")
    filtered_mt = mt.filter_cols(
        (mt.hgdp | mt.tgp) & ~hl.literal(populations).contains(mt.project_pop),
        keep=False,
    )
    filtered_mt = filtered_mt.filter_rows(hl.agg.any(filtered_mt.GT.is_non_ref()))
    logger.info(f"{filtered_mt.count_cols()} samples in filtered MatrixTable.")
    logger.info("Exporting sample map...")
    output_path = f"{output_path}/sample_map.tsv"
    all_ids = filtered_mt.cols()
    all_ids.export(output_path)

    return filtered_mt


def main(args):
    """Subset a matrix table to specified samples and across specified contigs."""
    contig = args.contig
    output_path = args.output_path

    if args.mt_path is not None:
        mt_path = args.mt_path
    elif args.mt_path is None and contig is not None:
        mt_path = MT_PATH.format(contig=contig)
    else:
        raise ValueError("Must provide either mt_path or contig.")
    mt = hl.read_matrix_table(mt_path)

    if args.test:
        logger.info("Filtering to test partitions...")
        mt = mt._filter_partitions(list(range(hl.eval(hl.min(mt.n_partitions(), 2)))))

    mt = subset_mt(mt, args.populations, output_path)

    output_path = (
        f"{output_path}/{contig}/{contig}_{args.subset_pop}.vcf.bgz"
        if contig
        else f"{output_path}/{args.subset_pop}.vcf.bgz"
    )

    hl.export_vcf(mt, output_path, tabix=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt-path", help="Path to the input MatrixTable.")
    parser.add_argument("--contig", help="Contig to subset.")
    parser.add_argument("--output-path", help="Path to the output bucket.")
    parser.add_argument(
        "--populations",
        nargs="+",
        help="List of populations to subset (e.g., ['afr', 'nfe']).",
        default=["afr", "nfe"],
    )
    parser.add_argument("--subset-pop", help="Population to subset.", default="afr")
    parser.add_argument(
        "--test",
        help="Subset to 2 partitions and variants on chr1.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing files.",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)


# testcase:
# mt_path = "gs://gnomad-lai/afr/sample_subsets/chr22/chr22_dense_biallelic_snps.mt"
# meta_path = "gs://gnomad-lai/afr/sample_subsets/afr_subset_samples_with_hgdp_tgp.ht"
# populations = ["afr", "nfe"]
# output_file_path = 'gs://kore-sandbox-storage/full_list.txt'
# filtered_mt = subsetPops(mt_path, meta_path, populations, output_file_path)


# complete the following steps
# 1) run subset_samples_and_variants function in gnomad_methods
# 2) export as vcf
# 3) run joint phasing with eagle

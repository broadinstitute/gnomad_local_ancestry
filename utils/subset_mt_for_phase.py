"""Module to subset MT for phasing."""

import logging
import hail as hl

from gnomad.utils.filtering import subset_samples_and_variants

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

MT_PATH = "gs://gnomad-lai/afr/sample_subsets/{contig}/{contig}_dense_biallelic_snps.mt"
META_PATH = "gs://gnomad-lai/afr/sample_subsets/afr_subset_samples_with_hgdp_tgp.ht"


def subset_mt(mt: hl.MatrixTable, meta: hl.Table, output_path: str):
    """
    Subset a MatrixTable based on specified populations and create a sample map.

    :param mt: Input MatrixTable.
    :param meta: Table of sample metadata.
    :param output_path: Path for the output sample map file.
    :return: Filtered MatrixTable.
    """
    logger.info("Annotating MatrixTable with meta information...")
    mt = mt.annotate_cols(**meta[mt.s])
    logger.info(f"{mt.count_cols()} samples in MatrixTable.")
    mt = subset_samples_and_variants(mt, f"{output_path}/full_sample_map.tsv")
    return mt


def get_subset_samples(populations: list, output_path: str) -> hl.Table:
    """
    Export samples to subset out of MTs.

    :param populations: List of populations to subset (e.g., ["afr", "nfe"]).
    :param output_path: Path for the output sample map file.
    :return: Table of samples to subset to.
    """
    meta = hl.read_table(META_PATH)
    meta = meta.select(
        hgdp=meta.subsets.hgdp,
        tgp=meta.subsets.tgp,
        project_pop=meta.project_meta.project_pop,
        project_subpop=meta.project_meta.project_subpop,
        hard_filters=meta.sample_filters.hard_filters,
    )

    logger.info("Removing HGDP/TGP samples from unspecified populations...")
    meta = meta.filter(
        ((meta.hgdp | meta.tgp) & ~hl.literal(populations).contains(meta.project_pop)),
        keep=False,
    )

    # Drop HGDP/TGP samples with hard fitlers (21 samples) as they were removed in the initial subsetting step.
    meta = meta.filter(
        (hl.is_defined(meta.hard_filters)) & (meta.hard_filters.length() > 0),
        keep=False,
    )

    meta = meta.annotate(
        sample_type=hl.if_else(
            (meta.hgdp | meta.tgp)
            & ((meta.project_subpop != "ASW") | (meta.project_subpop != "ACB")),
            "reference",
            "cohort",
        )
    )
    logger.info("Exporting full sample map...")
    meta.export(f"{output_path}/full_sample_map.tsv")

    logger.info("Exporting data type specific sample maps...")
    for sample_type in ["reference", "cohort"]:
        sample_map = meta.filter(meta.sample_type == sample_type)
        sample_map_path = f"{output_path}/{sample_type}_map.tsv"
        sample_map.export(sample_map_path)

    return meta


def main(args):
    """Subset a matrix table to specified samples and across specified contigs."""
    contigs = args.contigs
    contigs = [x for x in range(1, 23)] if args.contigs == ["all"] else contigs

    for contig in contigs:
        contig = "chr" + str(contig)
        output_path = args.output_path

        if args.mt_path is not None:
            mt_path = args.mt_path
        elif args.mt_path is None and contig is not None:
            mt_path = MT_PATH.format(contig=contig)
        else:
            raise ValueError("Must provide either mt_path or contig.")

        logger.info("Subsetting %s MatrixTable...", contig if contig else "passed")
        mt = hl.read_matrix_table(mt_path)
        logger.info(
            "MatrixTable loaded with %s variants across %s samples.",
            mt.count_rows(),
            mt.count_cols(),
        )

        if args.test:
            logger.info("Filtering to test partitions...")
            mt = mt._filter_partitions(
                list(range(hl.eval(hl.min(mt.n_partitions(), 2))))
            )
            output_path = output_path.replace("gnomad-lai", "gnomad-tmp-4day")

        logger.info("Getting samples to subset...")
        meta = get_subset_samples(args.populations, output_path)
        mt = subset_mt(mt, meta, output_path)

        output_path = f"{output_path}/{contig}" if contig else f"{output_path}"

        logger.info("Exporting MatrixTable and VCF...")
        mt = mt.checkpoint(
            f"{output_path}/{contig+'_' if contig else ''}{args.subset_pop}.mt",
            overwrite=args.overwrite,
        )

        hl.export_vcf(
            mt,
            f"{output_path}/{contig+'_' if contig else ''}{args.subset_pop}.vcf.bgz",
            tabix=True,
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mt-path", help="Path to the input MatrixTable.")
    parser.add_argument(
        "--contigs", help="Contig to subset.", nargs="+", default=["all"]
    )
    parser.add_argument(
        "--output-path",
        help="Path to the output bucket.",
        default="gs://gnomad-lai/afr/sample_subsets",
    )
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

"""Module to subset MT for phasing."""

import hail as hl


def subset_mt(mt_path: str, meta_path: str, populations: list, output_file_path: str):
    """
    Subset a MatrixTable based on specified populations and create a sample map.

    :param mt_path: Path to the input MatrixTable.
    :param meta_path: Path to the meta table.
    :param populations: List of populations to subset (e.g., ["afr", "nfe"]).
    :param output_file_path: Path for the output sample map file.
    :return: Filtered MatrixTable.
    """
    # Read MatrixTable and meta table
    mt = hl.read_matrix_table(mt_path)
    meta = hl.read_table(meta_path)

    # Annotate MatrixTable with meta information
    mt = mt.annotate_cols(**meta[mt.s])

    # Filter MatrixTable based on specified populations
    filtered_mt = mt.filter_cols(
        hl.if_else(
            hl.agg.any(mt.subsets.hgdp) | hl.agg.any(mt.subsets.tgp),
            hl.literal(populations).contains(mt.project_meta.project_pop),
            False,
        )
    )

    # Collect sample IDs
    ref_ids = list(filtered_mt.s.collect())
    individual_ids = list(mt.s.collect())

    # Combine ref_ids and individual_ids
    all_ids = set(ref_ids + individual_ids)

    # Write the sample map
    with hl.hadoop_open(output_file_path, "w") as f:
        # Write the header row with an 's'
        f.write("s\n")

        # Write each sample ID to the file
        for sample_id in all_ids:
            f.write(f"{sample_id}\n")

    return filtered_mt


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

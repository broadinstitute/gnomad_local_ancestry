# noqa: D100
import argparse
import hail as hl
import logging

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def import_lai_mt(
    anc: str,
    tractor_output_path: str = "tractor/test_path",
    file_extension: str = "",
    dosage: bool = True,
) -> hl.MatrixTable:
    """
    Import Tractor's dosage and hapcount files as hail MatrixTables.

    :param anc: File's ancestry.
    :param output_path: Path to Tractor's output files, defaults to "tractor/test_path".
    :param file_extension: If zipped, zip file extension, defaults to ''
    :param dosage: Whether the ancestry file being converted is dosage. When true, dosage file will be converted and when false the haps file will be converted. Defaults to True
    :return: hl.MatrixTable
    """
    row_fields = {
        "CHROM": hl.tstr,
        "POS": hl.tstr,
        "ID": hl.tstr,
        "REF": hl.tstr,
        "ALT": hl.tstr,
    }
    mt = hl.import_matrix_table(
        f"{tractor_output_path}.anc{anc}.{'dosage' if dosage else 'hapcount'}.txt{file_extension}",
        row_fields=row_fields,
        min_partitions=32,
    )
    mt = mt.key_rows_by().drop("row_id", "ID")

    variant = hl.parse_variant(
        mt.CHROM + ":" + mt.POS + ":" + mt.REF + ":" + mt.ALT, reference_genome="GRCh38"
    )
    mt = mt.key_rows_by(locus=variant["locus"], alleles=variant["alleles"]).drop(
        "CHROM", "POS", "REF", "ALT"
    )
    return mt


def generate_anc_mt_dict(
    ancs: dict, output_path: str = "tractor/test_path", file_extension: str = ""
) -> dict:
    """
    Generate dictionary where the key is ancestry and values are the ancestry's corresponding MatrixTable with hap and dosage annotations.

    :param ancs: Dictionary with keys as numerical value of msp file and values as the corresponding ancestry.
    :param output_path: Path to Tractor's output files, defaults to "tractor/test_path".
    :param file_extension: If zipped, zip file extension, defaults to ''.
    :return: Dictionary with ancestry(key) and corresponding Matrixtable(value).
    """
    logger.info(f"Generating ancestry matrixtable dictionary, ancestries are -> {ancs}")
    ancestry_mts = {}
    for num, anc in ancs.items():
        dos = import_lai_mt(num, output_path, file_extension, dosage=True)
        hap = import_lai_mt(num, output_path, file_extension, dosage=False)
        dos = dos.annotate_entries(
            **{f"{anc}_dos": dos.x, f"{anc}_hap": hap[dos.row_key, dos.col_id].x}
        )
        dos = dos.drop("x")
        ancestry_mts[anc] = dos
    return ancestry_mts


def get_msp_ancestries(msp_file: str = "tractor/test.msp.tsv") -> dict:
    """
    Parse msp header into dictionary of numeric keys and corresponding ancestry as values.

    :param msp_file: Path to msp file output by LAI tool like RFMixv2, defaults to "tractor/test.msp.tsv"
    :return: Dictionary of numeric keys and corresponding ancestry as values.
    """
    ancestries = ""
    with open(msp_file) as mspfile:
        line = mspfile.readline()
        if line.startswith("#Subpopulation order/codes:"):
            ancestries = {
                anc.split("=")[1]: anc.split("=")[0]
                for anc in line.strip().split(":")[1].strip().split("\t")
            }
    logger.info(f"Ancestries in msp file are {ancestries}")
    if len(ancestries) == 0:
        raise ValueError("Cannot find ancestries")
    return ancestries


def generate_joint_vcf(
    msp_file: str, tractor_output: str, output_path: str, is_zipped: bool = True
):
    """
    Generate a joint VCF from Trator's output files with ancestry-specific AC,AN,AF annotations.

    :param msp_file: Path to msp file output by LAI tool like RFMixv2, defaults to "tractor/test.msp.tsv"
    :param tractor_output_filepaths: Path to tractor output files without .hapcount.txt and .dosage.txt, e.g. /Tractor/output/test_run
    """
    logger.info(
        f"Generating joint VCF with annotated AFs. msp file is: {msp_file}, tractor output is {tractor_output}, is_zipped is {is_zipped}"
    )
    file_extension = ".gz" if is_zipped else ""
    ancestries = get_msp_ancestries(msp_file)
    anc_mts = generate_anc_mt_dict(
        ancs=ancestries, output_path=tractor_output, file_extension=file_extension,
    )
    entry_ancs = anc_mts.copy()
    anc, mt = entry_ancs.popitem()
    dos_hap_dict = {}
    callstat_dict = {}
    for anc, anc_mt in entry_ancs.items():
        dos_hap_dict.update(
            {
                f"{anc}_dos": anc_mt[mt.row_key, mt.col_key][f"{anc}_dos"],
                f"{anc}_hap": anc_mt[mt.row_key, mt.col_key][f"{anc}_hap"],
            }
        )

    mt = mt.annotate_entries(**dos_hap_dict)
    for anc, anc_mt in anc_mts.items():
        callstat_dict.update(
            {
                f"{anc}_AC": hl.agg.sum(mt[f"{anc}_dos"]),
                f"{anc}_AN": hl.agg.sum(mt[f"{anc}_hap"]),
                f"{anc}_AF": hl.if_else(
                    hl.is_defined(
                        hl.agg.sum(mt[f"{anc}_dos"]) / hl.agg.sum(mt[f"{anc}_hap"])
                    ),
                    hl.agg.sum(mt[f"{anc}_dos"]) / hl.agg.sum(mt[f"{anc}_hap"]),
                    0,
                ),
            }
        )
    mt = mt.annotate_rows(info=hl.struct(**callstat_dict))
    ht = mt.rows()
    hl.export_vcf(ht, f"{output_path}_lai_annotated.vcf.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--msp-file", help="Output from LAI program like RFMix_v2.", required=True
    )
    parser.add_argument(
        "--tractor-output",
        help="Path to tractor output files without .hapcount.txt and .dosage.txt, e.g. /Tractor/output/test_run",
        required=True,
    )
    parser.add_argument(
        "--output-path",
        help="Optional output path for files and file prefix, e.g. ~/test_data/test1 .",
    )
    parser.add_argument(
        "--is-zipped", help="Input files are gzipped.", action="store_true"
    )
    args = parser.parse_args()
    generate_joint_vcf(**vars(args))

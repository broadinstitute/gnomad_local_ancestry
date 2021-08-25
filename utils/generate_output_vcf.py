# noqa: D100
import argparse
import logging
from typing import Dict

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def import_lai_mt(
    anc: int,
    tractor_output_path: str = "tractor/test_path",
    file_extension: str = "",
    dosage: bool = True,
    min_partitions: int = 32,
) -> hl.MatrixTable:
    """
    Import Tractor's dosage and hapcount files as hail MatrixTables.

    :param anc: File's ancestry.
    :param output_path: Path to Tractor's output files, defaults to "tractor/test_path".
    :param file_extension: If zipped, zip file extension, defaults to ''.
    :param dosage: Whether the ancestry file being converted is a dosage file.
    When true, dosage file will be converted, and when false, the haps file will be converted. Defaults to True.
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :return: Dosage or hapcounts MatrixTable.
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
        min_partitions=min_partitions,
    )
    mt = mt.key_rows_by().drop("row_id", "ID")

    return mt.key_rows_by(
        locus=hl.locus(mt.CHROM, mt.POS), alleles=[mt.ref, mt.alt]
    ).drop("CHROM", "POS", "REF", "ALT")


def generate_anc_mt_dict(
    ancs: Dict[int, str],
    output_path: str = "tractor/test_path",
    file_extension: str = "",
    min_partitions: int = 32,
) -> Dict[str, hl.MatrixTable]:
    """
    Generate dictionary where the key is ancestry and values are the ancestry's corresponding MatrixTable with hap and dosage annotations.

    :param ancs: Dictionary with keys as numerical value of msp file and values as the corresponding ancestry.
    :param output_path: Path to Tractor's output files, defaults to "tractor/test_path".
    :param file_extension: If zipped, zip file extension, defaults to ''.
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :return: Dictionary with ancestry (key) and corresponding Matrixtable(value).
    """
    logger.info(
        "Generating ancestry matrixtable dictionary, ancestries are -> %s", ancs
    )
    ancestry_mts = {}
    for num, anc in ancs.items():
        dos = import_lai_mt(
            num, output_path, file_extension, dosage=True, min_partitions=min_partitions
        )
        hap = import_lai_mt(
            num,
            output_path,
            file_extension,
            dosage=False,
            min_partitions=min_partitions,
        )
        dos = dos.transmute_entries(
            **{f"{anc}_dos": dos.x, f"{anc}_hap": hap[dos.row_key, dos.col_id].x}
        )
        ancestry_mts[anc] = dos
    return ancestry_mts


def get_msp_ancestries(msp_file: str = "tractor/test.msp.tsv") -> Dict[int, str]:
    """
    Parse msp header into dictionary of numeric keys and corresponding ancestry strings as values.

    :param msp_file: Path to msp file output by LAI tool like RFMixv2, defaults to "tractor/test.msp.tsv".
    :return: Dictionary of numeric keys and corresponding ancestry as values.
    """
    ancestries = ""
    with open(msp_file) as mspfile:
        line = mspfile.readline()
    if line.startswith("#Subpopulation order/codes:"):
        # Header line of msp file is "#Subpopulation order/codes: ANC_I=0    ANC_J=1    ANC_K=2"
        ancestries = {
            anc.split("=")[1]: anc.split("=")[0]
            for anc in line.strip().split(":")[1].strip().split("\t")
        }
        logger.info("Ancestries in msp file are %s", ancestries)
        if len(ancestries) == 0:
            raise ValueError("Cannot find ancestries in header")
    return ancestries


def generate_joint_vcf(
    msp_file: str,
    tractor_output: str,
    output_path: str,
    is_zipped: bool = True,
    min_partitions: int = 32,
) -> None:
    """
    Generate a joint VCF from Trator's output files with ancestry-specific AC, AN, AF annotations.

    :param msp_file: Path to msp file output by LAI tool like RFMixv2, defaults to "tractor/test.msp.tsv".
    :param tractor_output_filepaths: Path to tractor output files without .hapcount.txt and .dosage.txt, e.g. /Tractor/output/test_run.
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :return: None; exports VCF to output path.
    """
    logger.info(
        "Generating joint VCF with annotated AFs. msp file is: %s, tractor output is %s, is_zipped is %s",
        msp_file,
        tractor_output,
        is_zipped,
    )
    file_extension = ".gz" if is_zipped else ""
    ancestries = get_msp_ancestries(msp_file)
    anc_mts = generate_anc_mt_dict(
        ancs=ancestries,
        output_path=tractor_output,
        file_extension=file_extension,
        min_partitions=min_partitions,
    )
    # Use one of the ancestry MTs as the base for the VCF export
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
    for anc in anc_mts:
        callstat_dict.update(
            {
                f"{anc}_AC": hl.agg.sum(mt[f"{anc}_dos"]),
                f"{anc}_AN": hl.agg.sum(mt[f"{anc}_hap"]),
                f"{anc}_AF": hl.if_else(
                    hl.agg.sum(mt[f"{anc}_hap"]) == 0,
                    0,
                    hl.agg.sum(mt[f"{anc}_dos"]) / hl.agg.sum(mt[f"{anc}_hap"]),
                ),
            }
        )
    ht = mt.annotate_rows(info=hl.struct(**callstat_dict)).rows()
    hl.export_vcf(ht, f"{output_path}_lai_annotated.vcf.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--msp-file", help="Output from LAI program like RFMix_v2.", required=True
    )
    parser.add_argument(
        "--tractor-output",
        help="Path to tractor output files without .hapcount.txt and .dosage.txt extensions, e.g. /Tractor/output/test_run",
        required=True,
    )
    parser.add_argument(
        "--output-path",
        help="Optional output path for files and file prefix, e.g. ~/test_data/test1 .",
    )
    parser.add_argument(
        "--is-zipped", help="Input files are gzipped.", action="store_true"
    )
    parser.add_argument(
        "--min-partitions",
        help="Minimum number of partitions to use when reading in tsv files as hail MTs, defaults to 32.",
        default=32,
    )
    args = parser.parse_args()
    generate_joint_vcf(**vars(args))

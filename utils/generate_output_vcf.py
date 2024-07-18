# noqa: D100
import argparse
import logging
from typing import Dict, List

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.filtering import filter_to_adj
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)

#command for converting file from vcf.gz to vcf.bgz
def bgzip_file(input_file: str, output_file: str) -> None:
    """
    Run bgzip on the input file to convert it to a block gzipped file.

    :param input_file: Path to the input gzip file.
    :param output_file: Path to the output bgzip file.
    """
    subprocess.run(f"gunzip -c {input_file} | bgzip -c > {output_file}", shell=True, check=True)
    
#convert file from vcf.gz to vcf.bgz
def convert_files_to_bgzip(tractor_output_path: str, ancestries: Dict[int, str], file_extension: str):
    for num in ancestries.keys():
        for file_type in ['dosage', 'hapcount']:
            input_file = f"{tractor_output_path}.anc{num}.{file_type}.txt{file_extension}"
            output_file = f"{tractor_output_path}.anc{num}.{file_type}.txt.bgz"
            logger.info(f"Converting {input_file} to {output_file}")
            bgzip_file(input_file, output_file)

def import_lai_mt(
    anc: int,
    tractor_output_path: str = "tractor/test_path",
    file_extension: str = "",
    dosage: bool = True,
    min_partitions: int = 32,
    batch_run: bool = True,
) -> hl.MatrixTable:
    """
    Import Tractor's dosage and hapcount files as hail MatrixTables.

    :param anc: File's ancestry.
    :param output_path: Path to Tractor's output files, defaults to "tractor/test_path".
    :param file_extension: If zipped, zip file extension, defaults to "".
    :param dosage: Whether the ancestry file being converted is a dosage file.
        When true, dosage file will be converted, and when false, haps file will be converted. Defaults to True.
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :param batch_run: Whether the run is for a batch run, defaults to True.
    :return: Dosage or hapcounts MatrixTable.
    """
    file_type = 'dosage' if dosage else 'hapcount'
    tractor_file = f"{tractor_output_path}.anc{anc}.{file_type}.txt{file_extension}"
    output_file = f"{tractor_output_path}.anc{anc}.{file_type}.table.ht"

    # return mt
    mt = hl.import_matrix_table(
        tractor_file,
        row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr},
        min_partitions=min_partitions,
        force_bgz=True
    )

    # Key rows by locus and alleles, and drop unnecessary fields
    mt = mt.key_rows_by(
        locus=hl.locus(mt.CHROM, mt.POS, reference_genome="GRCh38"),
        alleles=[mt.REF, mt.ALT]
    ).drop("CHROM", "POS", "REF", "ALT", "ID")
    sample_ids = mt.col_key.collect()
    mt = mt.rename({'col_id': 's'})
    mt = mt.drop('row_id')

    return mt


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
    :param file_extension: If zipped, zip file extension, defaults to "".
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :return: Dictionary with ancestry (key) and corresponding Matrixtable (value).
    """
    logger.info("Generating the ancestry matrixtable dictionary, ancestries are -> %s", ancs)
    ancestry_mts = {}
    for num, anc in ancs.items():
        dos = import_lai_mt(num, output_path, file_extension, dosage=True, min_partitions=min_partitions)
        hap = import_lai_mt(num, output_path, file_extension, dosage=False, min_partitions=min_partitions)

        # Transmute entries to include both dosage and hapcount
        dos = dos.transmute_entries(
            **{f"{anc}_dos": hl.int32(dos.x), f"{anc}_hap": hl.int32(hap[dos.row_key, dos.col_key].x)}
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
    mt_path_for_adj: str = "",
    add_gnomad_af: bool = False,
    gnomad_af_pops: List[str] = ["amr", "afr", "eas", "nfe"],
) -> None:
    """
    Generate a joint VCF from Trator's output files with ancestry-specific AC, AN, AF annotations.

    :param msp_file: Path to msp file output by LAI tool like RFMixv2, defaults to "tractor/test.msp.tsv".
    :param tractor_output_filepaths: Path to tractor output files without .hapcount.txt and .dosage.txt, e.g. /Tractor/output/test_run.
    :param min_partitions: Minimum partitions to use when reading in tsv files as hail MTs, defaults to 32.
    :param mt_path_for_adj: Path to MT to filter to high quality genotypes before calculating AC.
    :param add_gnomad_af: Add gnomAD's population AFs.
    param gnomad_af_pops: gnomAD continental pop's for AF annotation.
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
    convert_files_to_bgzip(tractor_output, ancestries, file_extension)
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

    if mt_path_for_adj:
        # This step requires access to an MT generated from pipeline's input VCF because Tractor's output does not contain the fields necessary for adj filtering (GT, GQ, DP, AB).
        logger.info("Filtering LAI output to adjusted genotypes...")
        adj_mt = hl.read_matrix_table(mt_path_for_adj)
        adj_mt = filter_to_adj(adj_mt)
        mt = mt.filter_entries(hl.is_defined(adj_mt[mt.row_key, mt.col_key]))

    for anc in anc_mts:
        logger.info("Calculating and annotating %s call stats", anc)
        callstat_dict.update(
            {
                f"AC_{anc}": hl.agg.sum(mt[f"{anc}_dos"]),
                f"AN_{anc}": hl.agg.sum(mt[f"{anc}_hap"]),
                f"AF_{anc}": hl.if_else(
                    hl.agg.sum(mt[f"{anc}_hap"]) == 0,
                    0,
                    hl.agg.sum(mt[f"{anc}_dos"]) / hl.agg.sum(mt[f"{anc}_hap"]),
                ),
            }
        )
    if add_gnomad_af:
        logger.info(
            "Annotating with gnomAD allele frequencies from %s pops...", gnomad_af_pops
        )
        gnomad_release = public_release("genomes").ht()
        callstat_dict.update(
            {
                f"gnomad_AF_{pop}": gnomad_release[mt.row_key].freq[
                    hl.eval(gnomad_release.freq_index_dict[f"{pop}-adj"])
                ]["AF"]
                for pop in gnomad_af_pops
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
        help="Path to tractor output files without .hapcount.txt and .dosage.txt extensions, e.g. /Tractor/output/test_run.",
        required=True,
    )
    parser.add_argument(
        "--output-path",
        help="Optional output path for files and file prefix, e.g. ~/test_data/test1.",
    )
    parser.add_argument(
        "--is-zipped", help="Input files are gzipped.", action="store_true"
    )
    parser.add_argument(
        "--min-partitions",
        help="Minimum number of partitions to use when reading in tsv files as hail MTs, defaults to 32.",
        default=32,
    )
    parser.add_argument(
        "--mt-path-for-adj",
        help="Path to hail MatrixTable generated from pipeline input VCF. Must contain GT, GQ, DP, and AB fields. If MT path provided, script will filter to high quality GTs only.",
    )

    parser.add_argument(
        "--add-gnomad-af",
        help="Add gnomAD population allele frequencies from AMR, NFE, AFR, and EAS.",
        action="store_true",
    )
    args = parser.parse_args()
    generate_joint_vcf(**vars(args))

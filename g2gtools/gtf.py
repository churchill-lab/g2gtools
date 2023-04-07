#
# Collection of functions related to GTF
#

# standard library imports
from bz2 import BZ2File
from collections import namedtuple
from collections import OrderedDict
from gzip import GzipFile
from typing import TextIO
import sys

# 3rd party library imports
# none

# local library imports
from .exceptions import G2GGTFError
from . import g2g
from . import g2g_utils
from . import vci

gtfInfoFields = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes",
]
GTFRecord = namedtuple("GTFRecord", gtfInfoFields)

ATTRIBUTES_TO_ALTER = {
    "gene_id": "gene_id",
    "transcript_id": "transcript_id",
    "exon_id": "exon_id",
    "protein_id": "protein_id",
    "ccds_id": "ccds_id",
}

ATTRIBUTES_TO_ALTER_GFF = {
    "ID": "ID",
    "Parent": "Parent",
    "Name": "Name",
    "Note": "Note",
}


class GTF(object):
    """
    Simple GTF object for parsing GTF files

    http://blogger.nextgenetics.net/?e=27

    Supports transparent gzip decompression.
    """

    def __init__(self, file_name: str):
        """
        Instantiate a GTF object.

        Args:
            file_name: The name of the GTF file.
        """
        self.file_name: str = file_name
        self.current_line: str = ""
        self.current_record: GTFRecord | None = None
        # self.reader: Iterator = iter(g2g_utils.open_resource(file_name))
        # self.reader = iter(g2g_utils.open_resource(file_name))
        self.reader: GzipFile | TextIO | BZ2File | None = g2g_utils.open_resource(
            file_name
        )

    def __iter__(self):
        return self

    def __next__(self):
        self.current_line = g2g_utils.s(self.reader.__next__())
        while self.current_line.startswith("#") or self.current_line.startswith("!"):
            self.current_line = g2g_utils.s(self.reader.__next__())

        self.current_record = parse_gtf_line(self.current_line)

        return self.current_record


def attributes_to_odict(attributes: str) -> OrderedDict[str, str]:
    """
    Parse the GTF attribute column and return a dictionary.

    Args:
        attributes: A string of attributes.

    Returns:
        A dictionary of attributes.
    """
    if attributes == ".":
        return OrderedDict()

    ret = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elem = attribute.strip().split(" ")
            key = elem[0]
            val = " ".join(elem[1:])

            if val[0] == '"':
                val = val[1:]
            if val[-1] == '"':
                val = val[0:-1]
            ret[key] = val

    return ret


def odict_to_attributes(attributes: dict[str, str]) -> str:
    """
    Parse the GTF attribute dictionary and return a string.

    Args:
        attributes: A dictionary of GTF attributes.

    Returns:
        A string representation.
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f'{k} "{v}"')
        temp_atts = "; ".join(atts)
        return temp_atts.rstrip() + ";"

    return "."


def parse_gtf_line(line: str) -> GTFRecord:
    """
    Parse the GTF/GFF line.

    Args:
        line: A line from GTF file.

    Returns:
        A GTFRecord.
    """
    elem = line.strip().split("\t")

    # If this fails, the file format is not standard-compatible
    if len(elem) != len(gtfInfoFields):
        raise G2GGTFError("Improperly formatted GTF file")

    data = {
        "seqid": None if elem[0] == "." else elem[0].replace('"', ""),
        "source": None if elem[1] == "." else elem[1].replace('"', ""),
        "type": None if elem[2] == "." else elem[2].replace('"', ""),
        "start": None if elem[3] == "." else int(elem[3]),
        "end": None if elem[4] == "." else int(elem[4]),
        "score": None if elem[5] == "." else float(elem[5]),
        "strand": None if elem[6] == "." else elem[6].replace('"', ""),
        "frame": None if elem[7] == "." else elem[7].replace('"', ""),
    }

    # since GTF will contain '"' vs. GFF will contain '=' in the elem[8]
    if '"' in elem[8]:
        data.update({"attributes": attributes_to_odict(elem[8])})
    elif "=" in elem[8]:
        data.update({"attributes": attributes_to_odict_gff(elem[8])})

    return GTFRecord(**data)


def convert_gtf_file(
    vci_file: str | vci.VCIFile,
    gtf_file_name_in: str,
    gtf_file_name_out: str | None = None,
    reverse: bool | None = False,
    debug_level: int | None = 0,
) -> None:
    """
    Convert GTF file.

    Args:
        vci_file: Name of the VCI file or a VCIFile object.
        gtf_file_name_in: Input GTF file to convert.
        gtf_file_name_out: Name of output GTF file, None for stdout.
        reverse: True to process VCI in reverse.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    logger = g2g.get_logger(debug_level)

    if isinstance(vci_file, vci.VCIFile):
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
        vci_file.parse(reverse)

    gtf_file_name_in = g2g_utils.check_file(gtf_file_name_in)
    logger.warn(f"GTF INPUT FILE: {gtf_file_name_in}")

    if gtf_file_name_out:
        output_file = g2g_utils.check_file(gtf_file_name_out, "w")
        unmapped_file = f"{output_file}.unmapped"
        gtf_out = open(output_file, "w")
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info("GTF OUTPUT FILE: {0}".format(output_file))
        logger.info("GTF UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(gtf_file_name_in)
        unmapped_file = f"{input_name}.unmapped"
        unmapped_file = g2g_utils.check_file(unmapped_file, "w")
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info("GTF OUTPUT: stdout")
        logger.info("GTF UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [""] if vci_file.is_haploid() else ["_L", "_R"]

    logger.info("Converting GTF file...")

    gtf_file = GTF(gtf_file_name_in)

    total = 0
    success = 0
    fail = 0

    # GTF is 1 based, bx-python is 0 based
    # when we do the querying, we subtract 1 from the GTF file start position
    # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
    #       bx-python (s, e) does it s <= x < e.

    for record in gtf_file:
        logger.debug(f"\nORIGINAL: {str(gtf_file.current_line).strip()}")

        total += 1

        if total % 100000 == 0:
            logger.info(f"Processed {total:,} lines")

        for lr in left_right:
            seq_id = f"{record.seqid}{lr}"
            mappings = vci_file.find_mappings(seq_id, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                logger.debug("Fail due to no mappings")
                gtf_unmapped_file.write(gtf_file.current_line)
                fail += 0
                continue
            else:
                logger.debug(f"{len(mappings)} mappings found")

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            logger.debug(f"({record.start - 1}, {record.end}) => ({start}, {end})")

            elem = gtf_file.current_line.rstrip().split("\t")
            elem[0] = seq_id
            elem[3] = start
            elem[4] = end

            if lr:
                attributes = attributes_to_odict(elem[8])
                for k, v in ATTRIBUTES_TO_ALTER.items():
                    if k in attributes:
                        attributes[k] = f"{attributes[k]}{lr}"
                elem[8] = odict_to_attributes(attributes)

            logger.debug("     NEW: {0}".format("\t".join(map(str, elem))))

            gtf_out.write("\t".join(map(str, elem)))
            gtf_out.write("\n")

    if gtf_out:
        gtf_out.close()

    if gtf_unmapped_file:
        gtf_unmapped_file.close()

    logger.warn(f"Generated {success:,} records from the original {total:,} records")
    logger.warn("GTF file converted")


def convert_gff_file(
    vci_file: str | vci.VCIFile,
    gff_file_name_in: str,
    gff_file_name_out: str | None = None,
    reverse: bool | None = False,
    debug_level: int | None = 0,
) -> None:
    """
    Convert GFF file.

    Args:
        vci_file: Name of the VCI file or a VCIFile object.
        gff_file_name_in: Input GFF file to convert.
        gff_file_name_out: Name of output GFF file, None for stdout.
        reverse: True to process VCI in reverse.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    logger = g2g.get_logger(debug_level)

    if isinstance(vci_file, vci.VCIFile):
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
        vci_file.parse(reverse)

    gff_file_name_in = g2g_utils.check_file(gff_file_name_in)
    logger.warn(f"GFF INPUT FILE: {gff_file_name_in}")

    if gff_file_name_out:
        gff_file_name_out = g2g_utils.check_file(gff_file_name_out, "w")
        unmapped_file = f"{gff_file_name_out}.unmapped"
        gff_out = open(gff_file_name_out, "w")
        gff_unmapped_file = open(unmapped_file, "w")
        logger.info(f"OUTPUT GFF FILE: {gff_file_name_out}")
        logger.info(f"UNMAPPED GFF FILE: {unmapped_file}")
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(gff_file_name_in)
        unmapped_file = f"{input_name}.unmapped"
        unmapped_file = g2g_utils.check_file(unmapped_file, "w")
        gff_out = sys.stdout
        gff_unmapped_file = open(unmapped_file, "w")
        logger.info("OUTPUT GFF: stdout")
        logger.info(f"UNMAPPED GFF FILE: {unmapped_file}")

    left_right = [""] if vci_file.is_haploid() else ["_L", "_R"]

    logger.info("Converting GFF file...")

    gff_file = GTF(gff_file_name_in)

    total = 0
    success = 0
    fail = 0

    # GTF is 1 based, bx-python is 0 based
    # when we do the querying, we subtract 1 from the GTF file start position
    # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
    #       bx-python (s, e) does it s <= x < e.

    for record in gff_file:
        logger.debug(f"ORIGINAL: {str(gff_file.current_line).strip()}")

        total += 1

        if total % 100000 == 0:
            logger.info(f"Processed {total:,} lines")

        for lr in left_right:
            seq_id = f"{record.seqid}{lr}"
            mappings = vci_file.find_mappings(seq_id, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                logger.debug("\tFail due to no mappings")
                gff_unmapped_file.write(gff_file.current_line)
                fail += 0
                continue
            else:
                logger.debug(f"{len(mappings)} mappings found")

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            logger.debug(f"({record.start}, {record.end}) => ({start}, {end})")

            elem = gff_file.current_line.rstrip().split("\t")
            elem[0] = seq_id
            elem[3] = start
            elem[4] = end

            if lr:
                attributes = attributes_to_odict_gff(elem[8])
                for k, v in ATTRIBUTES_TO_ALTER_GFF.items():
                    if k in attributes:
                        attributes[k] = f"{attributes[k]}{lr}"
                elem[8] = odict_to_attributes_gff(attributes)

            logger.debug("     NEW: {0}".format("\t".join(map(str, elem))))

            gff_out.write("\t".join(map(str, elem)))
            gff_out.write("\n")

    if gff_out:
        gff_out.close()

    if gff_unmapped_file:
        gff_unmapped_file.close()

    logger.warn(f"Generated {success:,} records from the original {total:,} records")
    logger.warn("GFF file converted")


def attributes_to_odict_gff(attributes: str) -> OrderedDict[str, str]:
    """
    Parse the GFF attribute column and return a dict.

    Args:
        attributes: The string of attributes.

    Returns:
        A dictionary of attributes.
    """
    if attributes == ".":
        return OrderedDict()

    ret_gff = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elem = attribute.strip().split("=")
            key = elem[0]
            val = ",".join(elem[1:])

            ret_gff[key] = val

    return ret_gff


def odict_to_attributes_gff(attributes: dict[str, str]) -> str:
    """
    Assemble the GFF attributes into a string.

    Args:
        attributes: A dictionary of attributes.

    Returns:
        A string version of the GFF attributes.
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f"{k}={v}")
        temp_atts = ";".join(atts)
        return temp_atts.rstrip()

    return "."

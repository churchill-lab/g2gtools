#
# Collection of functions related to GTF
#

# standard library imports
from collections import namedtuple
from collections import OrderedDict
import sys

# 3rd party library imports
# none

# local library imports
from .exceptions import G2GGTFError
from . import g2g
from . import g2g_utils
from . import vci

gtfInfoFields = [
    "seqid", "source", "type", "start", "end",
    "score", "strand", "frame", "attributes"
]
GTFRecord = namedtuple("GTFRecord", gtfInfoFields)

ATTRIBUTES_TO_ALTER = {
    "gene_id": "gene_id",
    "transcript_id": "transcript_id",
    "exon_id": "exon_id",
    "protein_id": "protein_id",
    "ccds_id": "ccds_id"
}

ATTRIBUTES_TO_ALTER_GFF = {
    "ID": "ID",
    "Parent": "Parent",
    "Name": "Name",
    "Note": "Note"
}


class GTF(object):
    """
    Simple GTF object for parsing GTF files

    http://blogger.nextgenetics.net/?e=27

    Supports transparent gzip decompression.
    """
    def __init__(self, file_name):
        if not file_name:
            raise G2GGTFError("A filename must be supplied")

        self.file_name = file_name
        self.current_line = None
        self.current_record = None
        self.reader = iter(g2g_utils.open_resource(file_name))

    def __iter__(self):
        return self

    def __next__(self):
        self.current_line = g2g_utils.s(self.reader.__next__())
        while self.current_line.startswith("#") or self.current_line.startswith("!"):
            self.current_line = g2g_utils.s(self.reader.__next__())

        self.current_record = parse_gtf_line(self.current_line)

        return self.current_record


def attributes_to_odict(attributes):
    """
    Parse the GTF attribute column and return a dict
    """
    if attributes == ".":
        return OrderedDict()

    ret = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elems = attribute.strip().split(" ")
            key = elems[0]
            val = " ".join(elems[1:])

            if val[0] == '"':
                val = val[1:]
            if val[-1] == '"':
                val = val[0:-1]
            ret[key] = val

    return ret


def odict_to_attributes(attributes):
    """
    Parse the GTF attribute column and return a dict
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f'{k} "{v}"')
        temp_atts = "; ".join(atts)
        return temp_atts.rstrip() + ";"

    return "."


def parse_gtf_line(line):  # works for both gtf and gff files
    """
    Parse the GTF line.

    :param line: a line from GTF file
    :type line: str
    :return:
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
        "frame": None if elem[7] == "." else elem[7].replace('"', '')}

    # since GTF will contain '"' vs. GFF will contain '=' in the elem[8]
    if '"' in elem[8]:
        data.update({
            "attributes": attributes_to_odict(elem[8])
        })
    elif "=" in elem[8]:
        data.update({
            "attributes": attributes_to_odict_gff(elem[8])
        })

    return GTFRecord(**data)


def convert_gtf_file(vci_file, input_file, output_file=None, reverse=False, debug_level=0):
    """
    Convert GTF coordinates.

    The mappings of coordinates are stored in the :class:`.chain.ChainFile` object.

    :param vci_file: VCI input file
    :type vci_file: :class:`.chain.ChainFile`
    :param input_file: the input GTF file
    :type input_file: string
    :param output_file: the output GTF file
    :type output_file: string
    :param reverse: reverse direction of original chain file
    :type reverse: boolean
    :return:
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

    input_file = g2g_utils.check_file(input_file)
    logger.warn(f"INPUT FILE: {input_file}")

    gtf_out = None
    gtf_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, "w")
        unmapped_file = f"{output_file}.unmapped"
        gtf_out = open(output_file, "w")
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info("OUTPUT FILE: {0}".format(output_file))
        logger.info("UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = f"{input_name}.unmapped"
        unmapped_file = g2g_utils.check_file(unmapped_file, "w")
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info("OUTPUT FILE: stdout")
        logger.info("UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [""] if vci_file.is_haploid() else ["_L", "_R"]

    logger.info("Converting GTF file...")

    gtf_file = GTF(input_file)

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
            seqid = f"{record.seqid}{lr}"
            mappings = vci_file.find_mappings(seqid, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                logger.debug("\tFail due to no mappings")
                gtf_unmapped_file.write(gtf_file.current_line)
                fail += 0
                continue
            else:
                logger.debug(f"{len(mappings)} mappings found")

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            logger.debug(f"({record.start-1}, {record.end}) => ({start}, {end})")

            elems = gtf_file.current_line.rstrip().split("\t")
            elems[0] = seqid
            elems[3] = start
            elems[4] = end

            if lr:
                attributes = attributes_to_odict(elems[8])
                for k,v in ATTRIBUTES_TO_ALTER.items():
                    if k in attributes:
                        attributes[k] = f"{attributes[k]}{lr}"
                elems[8] = odict_to_attributes(attributes)

            logger.debug("     NEW: {0}".format("\t".join(map(str, elems))))

            gtf_out.write("\t".join(map(str, elems)))
            gtf_out.write("\n")

    gtf_out.close()
    gtf_unmapped_file.close()

    logger.warn(f"Generated {success:,} records from the original {total:,} records")
    logger.warn("GTF file converted")


#####################################################################################################
# All the functions below this line are exclusively for GFF parsing.

def convert_gff_file(vci_file, input_file, output_file=None, reverse=False, debug_level=0):
    """
    Convert GFF coordinates.

    The mappings of coordinates are stored in the :class:`.chain.ChainFile` object.

    :param vci_file: VCI input file
    :type vci_file: :class:`.chain.ChainFile`
    :param input_file: the input GTF file
    :type input_file: string
    :param output_file: the output GTF file
    :type output_file: string
    :param reverse: reverse direction of original chain file
    :type reverse: boolean
    :return:
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

    input_file = g2g_utils.check_file(input_file)
    logger.warn(f"INPUT FILE: {input_file}")

    gtf_out = None
    gtf_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, "w")
        unmapped_file = f"{output_file}.unmapped"
        gtf_out = open(output_file, "w")
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info(f"OUTPUT FILE: {output_file}")
        logger.info(f"UNMAPPED FILE: {unmapped_file}")
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = f"{input_name}.unmapped"
        unmapped_file = g2g_utils.check_file(unmapped_file, "w")
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, "w")
        logger.info("OUTPUT FILE: stdout")
        logger.info(f"UNMAPPED FILE: {unmapped_file}")

    left_right = [""] if vci_file.is_haploid() else ["_L", "_R"]

    logger.info("Converting GFF file...")

    gff_file = GTF(input_file)

    total = 0
    success = 0
    fail = 0

    # GTF is 1 based, bx-python is 0 based
    # when we do the querying, we subtract 1 from the GTF file start position
    # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
    #       bx-python (s, e) does it s <= x < e.

    for record in gff_file:
        logger.debug(f"\nORIGINAL: {str(gff_file.current_line).strip()}")

        total += 1

        if total % 100000 == 0:
            logger.info(F"Processed {total:,} lines")

        for lr in left_right:
            seqid = f"{record.seqid}{lr}"
            mappings = vci_file.find_mappings(seqid, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                logger.debug("\tFail due to no mappings")
                gtf_unmapped_file.write(gff_file.current_line)
                fail += 0
                continue
            else:
                logger.debug(f"{len(mappings)} mappings found")

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            logger.debug(f"({record.start}, {record.end}) => ({start}, {end})")

            elems = gff_file.current_line.rstrip().split("\t")
            elems[0] = seqid
            elems[3] = start
            elems[4] = end

            if lr:
                attributes = attributes_to_odict_gff(elems[8])
                for k,v in ATTRIBUTES_TO_ALTER_GFF.items():
                    if k in attributes:
                        attributes[k] = f"{attributes[k]}{lr}"
                elems[8] = odict_to_attributes_gff(attributes)

            logger.debug("     NEW: {0}".format("\t".join(map(str, elems))))

            gtf_out.write("\t".join(map(str, elems)))
            gtf_out.write("\n")

    gtf_out.close()
    gtf_unmapped_file.close()

    logger.warn(f"Generated {success:,} records from the original {total:,} records")
    logger.warn("GFF file converted")


def attributes_to_odict_gff(attributes):
    """
    Parse the GFF attribute column and return a dict
    """

    if attributes == ".":
        return OrderedDict()

    ret_gff = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elems = attribute.strip().split("=")
            key = elems[0]
            val = ",".join(elems[1:])

            ret_gff[key] = val

    return ret_gff


def odict_to_attributes_gff(attributes):
    """
    Parse the GFF attribute column and return a dict
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f"{k}={v}")
        temp_atts = ";".join(atts)
        return temp_atts.rstrip()

    return "."


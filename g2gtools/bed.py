#
# Collection of functions related to BED files
#
# 0-based
#

# standard library imports
import collections
import sys

# 3rd party library imports
# none

# local library imports
from . import g2g
from . import g2g_utils
from .exceptions import G2GBedError
from .vci import VCIFile

bed_fields = ["chrom", "start", "end", "name", "score", "strand", "extra"]
BEDRecord = collections.namedtuple("BEDRecord", bed_fields)


class BED(object):
    """
    Simple BED object for parsing and iterating BED files.

    Supports transparent gzip decompression.
    """
    def __init__(self, filename):
        if not filename:
            raise G2GBedError("A filename must be supplied")

        self.filename = filename
        self.current_line = None
        self.current_line_is_bed = False
        self.current_record = None
        self.reader = g2g_utils.open_resource(filename)
        self.nitems = None
        self.current_line_no = 0

    def __iter__(self):
        return self

    def __next__(self):

        self.current_line = g2g_utils.s(self.reader.__next__())
        self.current_line_no += 1

        while self.current_line and len(self.current_line.strip()) == 0:
            self.current_line = self.reader.__next__()
            self.current_line_no += 1

        if self.current_line.startswith("track"):
            self.current_line = self.current_line.strip()
            self.current_line_is_bed = False
            self.current_record = None
            return None

        self.current_line_is_bed = True
        elem = self.current_line.strip().split("\t")

        if not self.nitems:
            self.nitems = len(elem)
        else:
            if self.nitems != len(elem):
                raise G2GBedError("Improperly formatted BED file")

        try:
            bed_data = {
                "chrom": elem[0],
                "start": int(elem[1]),
                "end": int(elem[2]),
                "name": elem[3] if self.nitems > 3 else None,
                "score": elem[4] if self.nitems > 4 else None,
                "strand": elem[5] if self.nitems > 5 else None,
                "extra": elem[6:] if self.nitems > 6 else None
            }

            self.current_record = BEDRecord(**bed_data)
            return self.current_record
        except IndexError as ie:
            raise G2GBedError((
                "Improperly formatted BED file, "
                f"line number: {self.current_line_no}, "
                f"line: {self.current_line}"
            ))
        except ValueError as ve:
            raise G2GBedError((
                "Improperly formatted BED file, "
                f"line number: {self.current_line_no}, "
                f"line: {self.current_line}"
            ))


# TODO: kb test
def convert_bed_file(vci_file, input_file, output_file=None, reverse=False, debug_level=0):
    """
    Convert BED coordinates.

    :param vci_file: VCI input file
    :type vci_file: :class:`.vci.VCIFile`
    :param input_file: the input BED file
    :type input_file: string
    :param output_file: the output BED file
    :type output_file: string
    :param reverse: reverse direction of original file
    :type reverse: boolean
    :return:
    """
    logger = g2g.get_logger(debug_level)

    if isinstance(vci_file, VCIFile):
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = VCIFile(vci_file)
        logger.warn(f"VCI FILE: {vci_file.filename}")
        logger.info(f"VCI FILE IS DIPLOID: {vci_file.is_diploid()}")
        vci_file.parse(reverse)

    input_file = g2g_utils.check_file(input_file)
    logger.warn(f"INPUT FILE: {input_file}")

    bed_out = None
    bed_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, "w")
        unmapped_file = f"{output_file}.unmapped"
        bed_out = open(output_file, "w")
        bed_unmapped_file = open(unmapped_file, "w")
        logger.info(f"OUTPUT FILE: {output_file}")
        logger.info(f"UNMAPPED FILE: {unmapped_file}")
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = f"{input_name}.unmapped"
        unmapped_file = g2g_utils.check_file(unmapped_file, "w")
        bed_out = sys.stdout
        bed_unmapped_file = open(unmapped_file, "w")
        logger.info("OUTPUT FILE: stdout")
        logger.info(f"UNMAPPED FILE: {unmapped_file}")

    left_right = [""] if vci_file.is_haploid() else ["_L", "_R"]

    logger.info("Converting BED file...")

    bed_file = BED(input_file)

    total = 0
    success = 0
    fail = 0

    for record in bed_file:

        logger.debug("\nORIGINAL: {0}".format(str(bed_file.current_line).strip()))

        total += 1

        if total % 100000 == 0:
            logger.info(f"Processed {total:,} lines")

        for lr in left_right:
            seqid = f"{record.chrom}{lr}"
            mappings = vci_file.find_mappings(
                seqid,
                record.start - 1,
                record.end
            )

            # unmapped
            if mappings is None:
                logger.debug("\tFail due to no mappings")
                bed_unmapped_file.write(bed_file.current_line)
                fail += 0
                continue
            else:
                logger.debug(f"{len(mappings)} mappings found")

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end
            elems = bed_file.current_line.rstrip().split("\t")

            logger.debug(f"({record.start-1}, {record.end}) => ({start}, {end})")
            logger.debug(elems)

            elems[0] = seqid
            elems[1] = start
            elems[2] = end

            temp_elems = "\t".join(map(str, elems))
            logger.debug(f"     NEW: {temp_elems}")

            bed_out.write(f"{temp_elems}\n")

    bed_out.close()
    bed_unmapped_file.close()

    logger.warn(f"Converted {success:,} of {total:,} records")
    logger.warn("BED file converted")


# -*- coding: utf-8 -*-

#
# Collection of functions related to BED files
#
# 0-based
#

from collections import namedtuple

from .chain import ChainFile
from .exceptions import G2GBedError, G2GLocationError
from .g2g_utils import get_logger
import g2g_fileutils as g2g_fu

LOG = get_logger()

bed_fields = ["chrom", "start", "end", "name", "score", "strand", "extra"]
BEDRecord = namedtuple("BEDRecord", bed_fields)


class BED(object):
    """
    Simple BED object for parsing BED files

    Supports transparent gzip decompression.
    """
    def __init__(self, filename):
        if not filename:
            raise G2GBedError("A filename must be supplied")

        self.filename = filename
        self.current_line = None
        self.current_line_is_bed = False
        self.current_record = None
        self.reader = g2g_fu.open_resource(filename)
        self.nitems = None
        self.current_line_no = 0

    def __iter__(self):
        return self

    def next(self):
        self.current_line = self.reader.next()
        self.current_line_no += 1

        while self.current_line and len(self.current_line.strip()) == 0:
            self.current_line = self.reader.next()
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
            bed_data = {'chrom': elem[0],
                        'start': int(elem[1]),
                        'end': int(elem[2]),
                        'name': elem[3] if self.nitems > 3 else None,
                        'score': elem[4] if self.nitems > 4 else None,
                        'strand':  elem[5] if self.nitems > 5 else None,
                        'extra': elem[6:] if self.nitems > 6 else None}

            self.current_record = BEDRecord(**bed_data)
            return self.current_record
        except IndexError, ie:
            LOG.debug(ie.message)
            raise G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))
        except ValueError, ve:
            LOG.debug(ve.message)
            raise G2GBedError("Improperly formatted BED file, line number: {0}, line: {1}".format(self.current_line_no, self.current_line))


def convert_bed_file(chain_file, input_file, output_file, reverse=False):
    """
    Convert BED coordinates.

    The mappings of coordinates are stored in the :class:`.chain.ChainFile` object.

    :param chain_file: chain file used for conversion
    :type chain_file: :class:`.chain.ChainFile`
    :param str file_in: the input BED file
    :type file_in: string
    :param file_out: the output BED file
    :type file_out: string
    :param reverse: reverse direction of original chain file
    :type reverse: boolean
    :return: Nothing
    """
    if not isinstance(chain_file, ChainFile):
        chain_file = g2g_fu.check_file(chain_file)

    input_file = g2g_fu.check_file(input_file)
    output_file_name = g2g_fu.check_file(output_file, 'w')
    unmapped_file_name = "{0}.unmapped".format(output_file_name)

    LOG.info("CHAIN FILE: {0}".format(chain_file))
    LOG.info("INPUT FILE: {0}".format(input_file))
    LOG.info("OUTPUT FILE: {0}".format(output_file_name))
    LOG.info("UNMAPPED FILE: {0}".format(unmapped_file_name))

    if not isinstance(chain_file, ChainFile):
        LOG.info("Parsing chain file...")
        chain_file = ChainFile(chain_file, reverse=reverse)
        LOG.info("Chain file parsed")

    bed_out = open(output_file_name, "w")
    bed_unmapped_file = open(unmapped_file_name, "w")

    LOG.info("Converting BED file")

    bed_file = BED(input_file)

    total = 0
    success = 0
    fail = 0

    # BED is 0 based, bx-python is 0 based

    try:
        for record in bed_file:
            # skip over "track" lines
            if not bed_file.current_line_is_bed:
                bed_out.write(bed_file.current_line)
                bed_out.write("\n")
                continue

            total += 1

            mappings = chain_file.find_mappings(record.chrom, record.start, record.end)

            # unmapped
            if mappings:
                success += 1
            else:
                LOG.debug("Fail due to no mappings")
                bed_unmapped_file.write(bed_file.current_line)
                fail += 1
                continue

            start = mappings[0].to_start
            end = mappings[-1].to_end

            LOG.debug("({0}, {1}) => ({2}, {3})".format(record.start, record.end, start, end))

            elems = bed_file.current_line.split()
            elems[1] = start
            elems[2] = end

            bed_out.write("\t".join(map(str, elems)))
            bed_out.write("\n")

        bed_out.close()
        bed_unmapped_file.close()

        LOG.info("Converted {0} of {1} records".format(success, total))
    except G2GLocationError, le:
        LOG.error("{0}: {1}".format(le.message, bed_file.current_line))

    LOG.info('BED file converted')


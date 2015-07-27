# -*- coding: utf-8 -*-

#
# Collection of functions related to GTF
#


from collections import namedtuple, OrderedDict
import os
import sys

from .g2g_utils import get_logger
from .exceptions import G2GGTFError, G2GValueError
from .chain import ChainFile

import g2g_fileutils as g2g_fu

gtfInfoFields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
GTFRecord = namedtuple('GTFRecord', gtfInfoFields)

LOG = get_logger()

class GTF(object):
    """
    Simple GTF object for parsing GTF files

    http://blog.nextgenetics.net/?e=27

    Supports transparent gzip decompression.
    """
    def __init__(self, file_name):
        if not file_name:
            raise G2GGTFError("A filename must be supplied")

        self.file_name = file_name
        self.current_line = None
        self.current_record = None
        self.reader = g2g_fu.open_resource(file_name)

    def __iter__(self):
        return self

    def next(self):
        self.current_line = self.reader.next()

        while self.current_line.startswith('#') or self.current_line.startswith('!'):
            self.current_line = self.reader.next()

        self.current_record = parse_gtf_line(self.current_line)

        return self.current_record


def parse_attributes(attributes):
    """
    Parse the GTF attribute column and return a dict
    """
    if attributes == ".":
        return OrderedDict()

    ret = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elems = attribute.strip().split(' ')
            key = elems[0]
            val = ' '.join(elems[1:])
            if val[0] == '"':
                val = val[1:]
            if val[-1] == '"':
                val = val[0:-1]
            ret[key] = val

    return ret


def parse_gtf_line(line):
    """
    Parse the GTF line.

    :param line: a line from GTF file
    :type line: str
    :return:
    """
    elem = line.strip().split("\t")

    # If this fails, the file format is not standard-compatible
    if len(elem) != len(gtfInfoFields):
        LOG.error(line)
        LOG.error("{0} != {1}".format(len(elem), len(gtfInfoFields)))
        raise G2GGTFError("Improperly formatted GTF file")

    data = {
        'seqid': None if elem[0] == '.' else elem[0].replace('"', ''),
        'source': None if elem[1] == '.' else elem[1].replace('"', ''),
        'type': None if elem[2] == '.' else elem[2].replace('"', ''),
        'start': None if elem[3] == '.' else int(elem[3]),
        'end': None if elem[4] == '.' else int(elem[4]),
        'score': None if elem[5] == '.' else float(elem[5]),
        'strand': None if elem[6] == '.' else elem[6].replace('"', ''),
        'frame': None if elem[7] == '.' else elem[7].replace('"', ''),
        'attributes': parse_attributes(elem[8])
    }

    return GTFRecord(**data)


def convert_gtf_file(chain_file, input_file, output_file, reverse=False):
    """
    Convert GTF coordinates.

    The mappings of coordinates are stored in the :class:`.chain.ChainFile` object.

    :param chain_file:
    :type chain_file: :class:`.chain.ChainFile`
    :param input_file: the input GTF file
    :type input_file: string
    :param output_file: the output GTF file
    :type output_file: string
    :param reverse: reverse direction of original chain file
    :type reverse: boolean
    :return:
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

    gtf_out = open(output_file_name, "w")
    gtf_unmapped_file = open(unmapped_file_name, "w")

    LOG.info("Converting GTF file...")

    gtf_file = GTF(input_file)

    total = 0
    success = 0
    fail = 0

    # GTF is 1 based, bx-python is 0 based
    # when we do the querying, we subtract 1 from the GTF file start position
    # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
    #       bx-python (s, e) does it s <= x < e.

    for record in gtf_file:

        LOG.debug("\nORIGINAL: {0}".format(str(gtf_file.current_line).strip()))

        total += 1

        if total % 100000 == 0:
            LOG.info("Processed {0:,} lines".format(total))

        mappings = chain_file.find_mappings(record.seqid, record.start - 1, record.end)

        # unmapped
        if mappings is None:
            LOG.debug("\tFail due to no mappings")
            gtf_unmapped_file.write(gtf_file.current_line)
            fail += 0
            continue
        else:
            LOG.debug("{0} mappings found".format(len(mappings)))

        success += 1
        start = mappings[0].to_start + 1
        end = mappings[-1].to_end

        LOG.debug("({0}, {1}) => ({2}, {3})".format(record.start - 1, record.end, start, end))

        elems = gtf_file.current_line.rstrip().split('\t')
        elems[3] = start
        elems[4] = end

        LOG.debug("     NEW: {0}".format("\t".join(map(str, elems))))

        gtf_out.write("\t".join(map(str, elems)))
        gtf_out.write("\n")

    gtf_out.close()
    gtf_unmapped_file.close()

    LOG.info("Converted {0:,} of {1:,} records".format(success, total))
    LOG.info('GTF file converted')


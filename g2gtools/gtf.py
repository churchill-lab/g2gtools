# -*- coding: utf-8 -*-

#
# Collection of functions related to GTF
#

from future.utils import implements_iterator

from collections import namedtuple, OrderedDict
import sys

from . import compat
from . import g2g
from . import g2g_utils
from . import vci
from . import exceptions

gtfInfoFields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
GTFRecord = namedtuple('GTFRecord', gtfInfoFields)

ATTRIBUTES_TO_ALTER = {'gene_id':'gene_id', 'transcript_id':'transcript_id', 'exon_id':'exon_id', 'protein_id':'protein_id', 'ccds_id':'ccds_id'}

ATTRIBUTES_TO_ALTER_GFF = {'ID':'ID', 'Parent':'Parent', 'Name':'Name', 'Note':'Note'}




LOG = g2g.get_logger()

@implements_iterator
class GTF(object):
    """
    Simple GTF object for parsing GTF files

    http://blog.nextgenetics.net/?e=27

    Supports transparent gzip decompression.
    """
    def __init__(self, file_name):
        if not file_name:
            raise exceptions.G2GGTFError("A filename must be supplied")

        self.file_name = file_name
        self.current_line = None
        self.current_record = None
        self.reader = iter(g2g_utils.open_resource(file_name))

    def __iter__(self):
        return self

    def __next__(self):
        if compat.is_py2:
            self.current_line = self.reader.next()
            while self.current_line.startswith('#') or self.current_line.startswith('!'):
                self.current_line = self.reader.next()
        else:
            self.current_line = g2g_utils.s(self.reader.__next__())
            while self.current_line.startswith('#') or self.current_line.startswith('!'):
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
            elems = attribute.strip().split(' ')
            key = elems[0]
            val = ' '.join(elems[1:])

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
            atts.append('{} "{}"'.format(k, v))
        temp_atts = "; ".join(atts)
        return temp_atts.rstrip() + ";"

    return '.'


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
        LOG.error(line)
        LOG.error("{0} != {1}".format(len(elem), len(gtfInfoFields)))
        raise exceptions.G2GGTFError("Improperly formatted GTF file")

    data = {
        'seqid': None if elem[0] == '.' else elem[0].replace('"', ''),
        'source': None if elem[1] == '.' else elem[1].replace('"', ''),
        'type': None if elem[2] == '.' else elem[2].replace('"', ''),
        'start': None if elem[3] == '.' else int(elem[3]),
        'end': None if elem[4] == '.' else int(elem[4]),
        'score': None if elem[5] == '.' else float(elem[5]),
        'strand': None if elem[6] == '.' else elem[6].replace('"', ''),
        'frame': None if elem[7] == '.' else elem[7].replace('"', '')}


    if '"' in elem[8]:  # since GTF will contain '"' vs. GFF will contain '=' in the elem[8]
        data.update({'attributes': attributes_to_odict(elem[8])})
    elif '=' in elem[8]: data.update({'attributes': attributes_to_odict_gff(elem[8])})

    return GTFRecord(**data)


def convert_gtf_file(vci_file, input_file, output_file=None, reverse=False):
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

    if isinstance(vci_file, vci.VCIFile):
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
        #vci_file.parse_fast(reverse)
        vci_file.parse(reverse)

        #LOG.info(vci_file.mapping_tree.keys())

    input_file = g2g_utils.check_file(input_file)
    LOG.info("INPUT FILE: {0}".format(input_file))

    gtf_out = None
    gtf_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, 'w')
        unmapped_file = "{0}.unmapped".format(output_file)
        gtf_out = open(output_file, "w")
        gtf_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: {0}".format(output_file))
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = "{0}.unmapped".format(input_name)
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: stdout")
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

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

        for lr in left_right:
            seqid = "{}{}".format(record.seqid, lr)
            mappings = vci_file.find_mappings(seqid, record.start - 1, record.end)

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
            elems[0] = seqid
            elems[3] = start
            elems[4] = end

            if lr:
                attributes = attributes_to_odict(elems[8])
                for k,v in ATTRIBUTES_TO_ALTER.items():
                    if k in attributes:
                        attributes[k] = '{}{}'.format(attributes[k], lr)
                elems[8] = odict_to_attributes(attributes)

            LOG.debug("     NEW: {0}".format("\t".join(map(str, elems))))

            gtf_out.write("\t".join(map(str, elems)))
            gtf_out.write("\n")

    gtf_out.close()
    gtf_unmapped_file.close()

    LOG.info("Generated {0:,} records from the original {1:,} records".format(success, total))
    LOG.info('GTF file converted')





#####################################################################################################
# All the functions below this line are exclusively for GFF parsing.

def convert_gff_file(vci_file, input_file, output_file=None, reverse=False):
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

    if isinstance(vci_file, vci.VCIFile):
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        LOG.info("VCI FILE: {0}".format(vci_file.filename))
        LOG.info("VCI FILE IS DIPLOID: {0}".format(vci_file.is_diploid()))
        #vci_file.parse_fast(reverse)
        vci_file.parse(reverse)

        #LOG.info(vci_file.mapping_tree.keys())

    input_file = g2g_utils.check_file(input_file)
    LOG.info("INPUT FILE: {0}".format(input_file))

    gtf_out = None
    gtf_unmapped_file = None

    if output_file:
        output_file = g2g_utils.check_file(output_file, 'w')
        unmapped_file = "{0}.unmapped".format(output_file)
        gtf_out = open(output_file, "w")
        gtf_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: {0}".format(output_file))
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(input_file)
        unmapped_file = "{0}.unmapped".format(input_name)
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, "w")
        LOG.info("OUTPUT FILE: stdout")
        LOG.info("UNMAPPED FILE: {0}".format(unmapped_file))

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']



    LOG.info("Converting GFF file...")

    gff_file = GTF(input_file)

    total = 0
    success = 0
    fail = 0

    # GTF is 1 based, bx-python is 0 based
    # when we do the querying, we subtract 1 from the GTF file start position
    # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
    #       bx-python (s, e) does it s <= x < e.

    for record in gff_file:
        LOG.debug("\nORIGINAL: {0}".format(str(gff_file.current_line).strip()))

        total += 1

        if total % 100000 == 0:
            LOG.info("Processed {0:,} lines".format(total))

        for lr in left_right:
            seqid = "{}{}".format(record.seqid, lr)
            mappings = vci_file.find_mappings(seqid, record.start - 1, record.end)

            # unmapped
            if mappings is None:
                LOG.debug("\tFail due to no mappings")
                gtf_unmapped_file.write(gff_file.current_line)
                fail += 0
                continue
            else:
                LOG.debug("{0} mappings found".format(len(mappings)))

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end


            LOG.debug("({0}, {1}) => ({2}, {3})".format(record.start - 1, record.end, start, end))

            elems = gff_file.current_line.rstrip().split('\t')
            elems[0] = seqid
            elems[3] = start
            elems[4] = end

            if lr:
                attributes = attributes_to_odict_gff(elems[8])
                for k,v in ATTRIBUTES_TO_ALTER_GFF.items():
                    if k in attributes:
                        attributes[k] = '{}{}'.format(attributes[k], lr)
                elems[8] = odict_to_attributes_gff(attributes)

            LOG.debug("     NEW: {0}".format("\t".join(map(str, elems))))

            gtf_out.write("\t".join(map(str, elems)))
            gtf_out.write("\n")

    gtf_out.close()
    gtf_unmapped_file.close()

    LOG.info("Generated {0:,} records from the original {1:,} records".format(success, total))
    LOG.info('GFF file converted')


def attributes_to_odict_gff(attributes):
    """
    Parse the GFF attribute column and return a dict
    """

    if attributes == ".":
        return OrderedDict()

    ret_gff = OrderedDict()
    for attribute in attributes.strip().split(";"):
        if len(attribute):
            elems = attribute.strip().split('=')
            key = elems[0]
            val = ','.join(elems[1:])

            ret_gff[key] = val

    return ret_gff


def odict_to_attributes_gff(attributes):
    """
    Parse the GFF attribute column and return a dict
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append('{}={}'.format(k, v))
        temp_atts = ";".join(atts)
        return temp_atts.rstrip()

    return '.'


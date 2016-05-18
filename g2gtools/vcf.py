# -*- coding: utf-8 -*-

#
# Collection of functions related to VCF files
#
# 1 based


from collections import namedtuple

from .g2g_utils import get_logger
from .g2g_fileutils import open_resource
from .exceptions import G2GVCFError

VCF_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLES']
VCFRecord = namedtuple('VCFRecord', VCF_FIELDS)

GT_DATA_FIELDS = ['ref', 'left', 'right', 'gt', 'fi', 'phase']
GTData = namedtuple('GTData', GT_DATA_FIELDS)

GENOTYPE_UNPHASED = '/'
GENOTYPE_PHASED = '|'

LOG = get_logger()

class VCF(object):
    """
    Simple VCF object for parsing VCF files
    """
    def __init__(self, file_name):
        if not file_name:
            raise G2GVCFError("A filename must be supplied")

        self.file_name = file_name
        self.samples = None
        self.current_line = None
        self.current_record = None
        self.reader = open_resource(file_name)
        self._parse_header()

    def _parse_header(self):
        self.current_line = self.reader.next()

        while self.current_line.startswith('##'):
            self.current_line = self.reader.next()

        if self.current_line.startswith('#'):
            elems = self.current_line.strip().split('\t')
            samples = elems[9:]
            self.samples = dict(zip(samples, (x for x in xrange(len(samples)))))
        else:
            raise G2GVCFError("Improperly formatted VCF file")

    def parse_gt(self, sample):
        if not sample:
            raise G2GVCFError("Sample must contain a value")

        sample_index = self.get_sample_index(sample)
        return parse_gt(self.current_record, sample_index)

    def __iter__(self):
        return self

    def next(self):
        self.current_line = self.reader.next()

        while self.current_line.startswith("#"):
            self.current_line = self.reader.next()

        self.current_record = parse_vcf_line(self.current_line)

        return self.current_record

    def get_sample_index(self, sample):
        if not sample:
            raise G2GVCFError("Sample must contain a value")

        if sample in self.samples:
            return self.samples[sample]

        parse_vcf_line()

        raise G2GVCFError("Unknown sample: '{0}'".format(sample))


def parse_vcf_line(line):
    """
    Parse a line in the VCF file.

    :param line: a line from the VCF file
    :type line: str
    :return: :class:`.vcf.VCFRecord`
    """

    if isinstance(line, basestring):
        if line.startswith('#'):
            return None

        elem = line.strip().split('\t')
    elif isinstance(line, list):
        elem = line

    try:
        quality = int(elem[5])
    except ValueError:
        try:
            quality = float(elem[5])
        except ValueError:
            quality = None

    filter_field = None

    if elem[6] != '.':
        filter_field = elem[6].split(';')

    info = elem[7]

    try:
        fmt = elem[8]
    except IndexError:
        fmt = None
    else:
        if fmt == '.':
            fmt = None

    return VCFRecord(elem[0], int(elem[1]), None if elem[2] == '.' else elem[2], elem[3],
                     elem[4].split(','), quality, filter_field, info, fmt, elem[9:])


def parse_gt(vcf_record, sample_index):
    """
    Parse the GT field within the VCF line.

    :param vcf_record: the VCF record
    :type vcf_record: :class:`.vcf.VCFRecord`
    :param sample_index: the strain or sample index
    :type sample_index: int
    :return: :class:`.vcf.GTData`
    """
    if sample_index is None:
        raise G2GVCFError("Sample index must contain a value")

    sample_data = vcf_record.SAMPLES[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None

    if sample_data != '.':
        gt_index = vcf_record.FORMAT.split(':').index('GT')
        fi_index = vcf_record.FORMAT.split(':').index('FI')

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = map(int, gt.split(GENOTYPE_PHASED))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = map(int, gt.split(GENOTYPE_UNPHASED))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError("Unknown phase in GT, {0}".format(gt))

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.REF
                else:
                    left = vcf_record.ALT[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_record.REF
                else:
                    right = vcf_record.ALT[genotypes[1]-1]

        except ValueError, ve:
            LOG.debug(ve)
        except IndexError, ie:
            LOG.debug(ie)
        try:
            fi = sample_data.split(':')[fi_index]
        except ValueError, ve:
            LOG.debug(ve)
        except IndexError, ie:
            LOG.debug(ie)

    return GTData(vcf_record.REF, left, right, gt, fi, phase)


def parse_gt_new(vcf_tuple, sample_index):
    """
    Parse the GT field within the VCF line.

    :param vcf_record: the VCF record
    :type vcf_record: :class:`.vcf.VCFRecord`
    :param sample_index: the strain or sample index
    :type sample_index: int
    :return: :class:`.vcf.GTData`
    """
    if sample_index is None:
        raise G2GVCFError("Sample index must contain a value")

    sample_data = vcf_tuple[sample_index]

    gt = None
    fi = None
    left = None
    right = None
    phase = None

    if sample_data != '.':
        gt_index = vcf_tuple.format.split(':').index('GT')

        try:
            fi_index = vcf_tuple.format.split(':').index('FI')
        except:
            fi_index = None

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = map(int, gt.split(GENOTYPE_PHASED))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = map(int, gt.split(GENOTYPE_UNPHASED))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError("Unknown phase in GT, {0}".format(gt))

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_tuple.ref
                else:
                    left = vcf_tuple.alt.split(',')[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_tuple.ref
                else:
                    right = vcf_tuple.alt.split(',')[genotypes[1]-1]

        except ValueError, ve:
            LOG.debug(ve)
        except IndexError, ie:
            LOG.debug(ie)

        if fi_index is not None:
            try:
                fi = sample_data.split(':')[fi_index]
            except ValueError, ve:
                LOG.debug(ve)
            except IndexError, ie:
                LOG.debug(ie)
        else:
            fi = None

    return GTData(vcf_tuple.ref, left, right, gt, fi, phase)

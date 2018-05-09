# -*- coding: utf-8 -*-

#
# Collection of functions related to VCF files
#
# 1 based

from future.utils import lmap
from past.builtins import xrange

from collections import namedtuple

import re

from . import g2g
from . import g2g_utils
from . import exceptions

VCF_FIELDS = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'samples']
VCFRecord = namedtuple('VCFRecord', VCF_FIELDS)

GT_DATA_FIELDS = ['ref', 'left', 'right', 'gt', 'fi', 'phase', 'gt_left', 'gt_right', 'is_snp']
GTData = namedtuple('GTData', GT_DATA_FIELDS)

GENOTYPE_UNPHASED = '/'
GENOTYPE_PHASED = '|'

REGEX_ALT = re.compile("(^[A|C|G|T]+)")

LOG = g2g.get_logger()


class VCFFile(object):
    """
    Simple VCF object for parsing VCF files
    """
    def __init__(self, file_name):
        if not file_name:
            raise exceptions.G2GVCFError("A filename must be supplied")

        self.file_name = file_name
        self.samples = None
        self.current_line = None
        self.current_record = None
        self.reader = g2g_utils.open_resource(file_name)
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
            raise exceptions.G2GVCFError("Improperly formatted VCF file")

    def parse_gt(self, sample):
        if sample is None:
            raise exceptions.G2GVCFError("Sample must contain a value")

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
        if sample is None:
            raise exceptions.G2GVCFError("Sample must contain a value")

        if sample in self.samples:
            return self.samples[sample]

        parse_vcf_line()

        raise exceptions.G2GVCFError("Unknown sample: '{0}'".format(sample))


def parse_vcf_line(line):
    """
    Parse a line in the VCF file.

    :param line: a line from the VCF file
    :type line: str
    :return: :class:`.vcf.VCFRecord`
    """

    if isinstance(line, str):
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
        raise exceptions.G2GVCFError("Sample index must contain a value")

    sample_data = vcf_record.samples[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find('<') == -1 and sample_data != '.':
    #if sample_data != '.':
        gt_index = vcf_record.format.split(':').index('GT')
        fi_index = vcf_record.format.split(':').index('FI')

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = lmap(int, gt.split(GENOTYPE_PHASED))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = lmap(int, gt.split(GENOTYPE_UNPHASED))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError("Unknown phase in GT, {0}".format(gt))

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt[genotypes[1]-1]

                gt_left = genotypes[0]
                gt_right = genotypes[1]

                # check for to see if ALT is <CN*> or something not ACGT
                if not REGEX_ALT.match(gt_left):
                    left = None
                    gt_left = None

                if not REGEX_ALT.match(gt_right):
                    right = None
                    gt_right = None

        except ValueError as ve:
            LOG.debug(ve)
        except IndexError as ie:
            LOG.debug(ie)
        try:
            fi = sample_data.split(':')[fi_index]
        except ValueError as ve:
            LOG.debug(ve)
        except IndexError as ie:
            LOG.debug(ie)

    is_snp = len(vcf_record.REF) == 1 == (len(left) if left else 0) == (len(right) if right else 0)
    return GTData(vcf_record.REF, left, right, gt, fi, phase, gt_left, gt_right, is_snp)


def parse_gt_tuple(vcf_record, sample_index):
    """
    Parse the GT field within the VCF line.
    """
    if sample_index is None:
        raise exceptions.G2GVCFError("Sample index must contain a value")

    sample_data = vcf_record[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find('<') == -1 and sample_data != '.':
        formats = vcf_record.format.split(':')
        gt_index = formats.index('GT')
        fi_index = formats.index('FI') if 'FI' in formats else None

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = lmap(int, gt.split(GENOTYPE_PHASED))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = lmap(int, gt.split(GENOTYPE_UNPHASED))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError("Unknown phase in GT, {0}".format(gt))

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt.split(',')[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt.split(',')[genotypes[1]-1]

                gt_left = genotypes[0]
                gt_right = genotypes[1]

                # check for to see if ALT is <CN*> or something not ACGT
                #if not REGEX_ALT.match(left) or not REGEX_ALT.match(right):
                #    LOG.error("VFC2VCI CN FOUND")
                #    gt = None
                #    fi = None
                #    left = None
                #    right = None
                #    phase = None
                #    gt_left = None
                #    gt_right = None

        except ValueError as ve:
            LOG.debug(ve)
        except IndexError as ie:
            LOG.debug(ie)
        try:
            if fi_index:
                fi = sample_data.split(':')[fi_index]
        except ValueError as ve:
            LOG.debug(ve)
        except IndexError as ie:
            LOG.debug(ie)

    is_snp = len(vcf_record.ref) == 1 == (len(left) if left else 0) == (len(right) if right else 0)
    return GTData(vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp)

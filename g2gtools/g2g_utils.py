# -*- coding: utf-8 -*-

#
# Collection of g2g utility functions and classes
#

import logging
import os
import re
import sys

from itertools import izip_longest
from string import maketrans

from . import G2GLocationError, G2GValueError

REGEX_LOCATION = re.compile("(\w*)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
REGEX_LOCATION_CHR = re.compile("(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
BASES = re.compile(r"([ACTGNactgn]+)")
TRANS = maketrans('ACTGNactgn', 'TGACNtgacn')

LOG = logging.getLogger('G2G')


class Location(object):
    def __init__(self, seqid, start=None, end=None, strand='+', name=None, original_base=0):
        if start:
            end = end if end else start + 1
        self.seqid = seqid
        self._strand = None
        self._start = None
        self._end = None
        self.original_base = original_base
        self.strand = strand
        self.start = start
        self.end = end
        self.name = name

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value):
        if value is not None:
            if value in ('+', '1', 1):
                self._strand = '+'
            elif value in ('-', '-1', -1):
                self._strand = '-'
            else:
                raise G2GLocationError("Illegal value for strand {0}, must be +, -, 1, -1".format(value))

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if value is not None:
            if value < 0:
                raise G2GLocationError("Illegal value for start {0}, start must be >= 0".format(value))
            self._start = value

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if value is not None:
            if self.start and value <= self.start:
                raise G2GLocationError("Illegal value for end {0}, end must be >= start {1}".format(value, self.start))

            self._end = value
        else:
            if self.start:
                self._end = self.start + 1

    def __str__(self):
        if self.original_base == 1:
            return "{0}:{1}-{2} ({3})".format(self.seqid, self.start+1, self.end, self.strand)
        else:
            return "{0}:{1}-{2} ({3})".format(self.seqid, self.start, self.end, self.strand)

    def __repr__(self):
        return "{0}:{1}-{2} ({3})".format(self.seqid, self.start, self.end, self.strand)


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


class G2GFormatter(logging.Formatter):

    err_fmt = "ERROR: %(msg)s"
    #dbg_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    dbg_fmt = "DEBUG: %(asctime)s: %(msg)s"
    info_fmt = "%(msg)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = G2GFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = G2GFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = G2GFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


def get_logger():
    return LOG


def configure_logging(level):
    """

    :param level: 0=ERROR, 1=INFO, 2+=DEBUG
    """
    global LOG
    LOG = logging.getLogger('G2G')

    handler = logging.StreamHandler()

    mfmt = G2GFormatter()
    handler.setFormatter(mfmt)

    if level == 0:
    #    LOG.setLevel(logging.ERROR)
    #elif level == 1:
        LOG.setLevel(logging.INFO)
    else:
        LOG.setLevel(logging.DEBUG)

    LOG.addHandler(handler)


def app_exit(message='', parser=None):
    """
    Print message, help, and exit.

    :param message: The message to print
    :param parser: the argument parser
    """
    if parser:
        parser.error(message)
    else:
        if message:
            print str(message)

    sys.exit()


def format_time(start, end):
    """
    Format length of time between start and end.

    :param start: the start time
    :param end: the end time
    :return: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def intersect_regions(chr1, start1, end1, chr2, start2, end2):
    """
    Return the intersection of two regions

    :param chr1: chromosome 1
    :param start1: start position 1
    :param end1: end position 1
    :param chr2: chromosome 2
    :param start2: start position 2
    :param end2: end position 2
    :return: the intersection of the coordinated
    """

    if chr1 != chr2:
        return None

    if int(start1) > int(end2) or int(end1) < int(start2):
        return None

    if int(start1) > int(end1) or int(start2) > int(end2):
        raise G2GLocationError("Start cannot be larger than end")

    return chr1, max(start1, start2), min(end1, end2)


def parse_chromosome(chrom_file):
    '''
    Parse the chromosome for length
    '''
    chromosomes = {}

    try:
        if not os.path.exists(chrom_file):
            raise IOError(chrom_file)

        fd = open(chrom_file, 'r')
        for line in fd:
            elem = line.strip().split()
            try:
                chromosome = elem[0]
                chromosome_length = int(elem[1])
                chromosomes[chromosome] = chromosome_length
            except ValueError, e1:
                pass # most likely header line
        fd.close()
    except IOError, e:
        message = "Error parsing chromosome files: {0}".format(e.message)
        print message
        return {}

    return chromosomes


def get_multiplier(factor):
    """
    Convert the factor into a number.

    :param factor: the string 'mb', 'm', or 'k'
    :return: 10000000, 1000000, 1000 or 1
    """
    if factor.lower() == 'mb':
        return 10000000
    elif factor.lower() == 'm':
        return 1000000
    elif factor.lower() == 'k':
        return 1000

    return 1


def parse_location(location_str, base=0):
    """
    Parse a string and return a location.

    Format expected is seqid|start|end, where | is some sort of delimiter.

    If end is None, 1 based was requested.

    :param location_str: a string representing a genomic location
    :return: a Location
    """
    loc_match = REGEX_LOCATION.match(location_str.strip().replace(',', ''))
    one_base = False

    if loc_match:
        loc_groups = loc_match.groups()
        identifier = loc_groups[0]
        start_base = loc_groups[2]
        start_mult = loc_groups[3]
        end_base = loc_groups[5]
        end_mult = loc_groups[6]

        if not identifier:
            raise G2GLocationError("Cannot parse location")

        if start_base:
            try:
                start_base = int(start_base)
            except ValueError:
                raise G2GLocationError("Start position is not numeric")
        else:
            raise G2GLocationError("Cannot parse start position")

        start = start_base

        if start_mult:
            if start_mult.lower() not in ['mb', 'm', 'k']:
                raise G2GLocationError("Unknown quantifier '{0}".format(start_mult))
            start = start_base * get_multiplier(start_mult)
        else:
            if base == 1:
                start -= 1

        if end_base:
            try:
                end_base = int(end_base)
            except ValueError:
                raise G2GLocationError("End position is not numeric")
        else:
            # special condition to allow for only 1 location to be retrieved
            # or from this start base on

            if not loc_groups[4] and not end_base and not end_mult:
                # only 1 position requested
                one_base = True
            else:
                raise G2GLocationError("Cannot parse end position")

        if end_mult:
            if end_mult.lower() not in ['mb', 'm', 'k']:
                raise G2GLocationError("Unknown quantifier '{0}".format(end_mult))
    else:
        raise G2GLocationError("Cannot parse location '{0}'".format(location_str))

    if one_base:
        end = None
    else:
        end = end_base * get_multiplier(end_mult)

    LOG.debug("identifier={0}, start={1}, end={2}, base={3}".format(identifier, start, end, base))

    return Location(identifier, start, end, '+', original_base=base)


def wrap_sequence(sequence, n=60, fillvalue=''):
    args = [iter(sequence)] * n
    for line in izip_longest(fillvalue=fillvalue, *args):
        yield ''.join(line + ('\n',))


def reverse_sequence(sequence):
    """
    Return the reverse of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse string or Sequence object
    """
    if isinstance(sequence, str):
        return sequence[::-1]
    else:
        raise ValueError("Cannot complement object of {0}, expecting string or Sequence object".format(type(sequence)))


def complement_sequence(sequence):
    """
    Return the complement of sequence

    :param sequence: either a string or Sequence object
    :return: the complement string or Sequence object
    """
    val = str(sequence)

    if len(re.findall(BASES, val)) == 1:  # finds invalid characters if > 1
        val = str(sequence).translate(TRANS)
    else:
        matches = re.findall(BASES, val)
        position = len(matches[0])
        raise ValueError("Sequence contains non-DNA character '{0}' at position {1:n}\n".format(val[position], position + 1))

    if isinstance(sequence, str):
        return val
    else:
        raise ValueError("Cannot complement object of {0}, expecting string or Sequence object".format(type(sequence)))


def reverse_complement_sequence(sequence):
    """
    Return the reverse-complement of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse-complement string or Sequence object
    """
    return reverse_sequence(complement_sequence(sequence))


if __name__ == '__main__':
    l = parse_location("1:0-33")
    print l
    print l.__repr__()

    l = parse_location("1:1-33")
    print l
    print l.__repr__()

    print '1'
    l = Location(1, 22, None, '-')
    print l
    print l.__repr__()

    print 'tryng to set l.start'
    l.start = 2
    print l
    print l.__repr__()

# -*- coding: utf-8 -*-

from io import StringIO

import logging
import re
import sys

from . import exceptions

REGEX_REGION = re.compile("(\w*)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
REGEX_REGION_CHR = re.compile("(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)

LOG = logging.getLogger('G2G')

#
#
# Logging
#
#


class G2GFormatter(logging.Formatter):

    err_fmt = "[g2gtools] %(msg)s"
    # err_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    # dbg_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    dbg_fmt = "[g2gtools debug] %(msg)s"
    info_fmt = "[g2gtools] %(msg)s"

    def __init__(self, fmt="[g2gtools] %(msg)s"):
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
    """
    Get the logger
    :return: logger object
    """
    return LOG


def configure_logging(level):
    """

    :param level: 0=ERROR, 1=INFO, 2+=DEBUG
    """
    global LOG
    LOG = logging.getLogger('G2G')

    handler = logging.StreamHandler(sys.stderr)

    mfmt = G2GFormatter()
    handler.setFormatter(mfmt)

    if level == 0:
        LOG.setLevel(logging.INFO)
    else:
        LOG.setLevel(logging.DEBUG)

    LOG.addHandler(handler)


def exit(message='', parser=None):
    """
    Print message, help, and exit.

    :param message: The message to print
    :param parser: the argument parser
    """
    if parser:
        parser.error(str(message))
    else:
        if message:
            sys.stderr.write("[g2gtools] Error: {}\n".format(str(message)))

    sys.exit(1)

#
#
# Region
#
#


class Region(object):
    """

    """
    def __init__(self, seq_id, start=None, end=None, strand='+', name=None, original_base=0):
        if start:
            end = end if end else start + 1
        self.seq_id = seq_id
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
                raise exceptions.G2GRegionError("Illegal value for strand {0}, must be +, -, 1, -1".format(value))

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if value is not None:
            if value < 0:
                raise exceptions.G2GRegionError("Illegal value for start {0}, start must be >= 0".format(value))
            self._start = value

    def get_start(self):
        if self.original_base == 1:
            return self.start
        return self.start - 1

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if value is not None:
            if self.start and value <= self.start:
                raise exceptions.G2GRegionError("Illegal value for end {0}, end must be >= start {1}".format(value, self.start))

            self._end = value
        else:
            if self.start:
                self._end = self.start + 1

    def str(self):
        if self.original_base == 1:
            if self.name:
                return "{} {}:{}-{} ({})".format(self.name, self.seq_id, self.start+1, self.end, self.strand)
            else:
                return "{}:{}-{} ({})".format(self.seq_id, self.start+1, self.end, self.strand)
        else:
            if self.name:
                return "{} {}:{}-{} ({})".format(self.name, self.seq_id, self.start, self.end, self.strand)
            else:
                return "{}:{}-{} ({})".format(self.seq_id, self.start, self.end, self.strand)

    def __str__(self):
        if self.original_base == 1:
            if self.name:
                return "{} {}:{}-{}".format(self.name, self.seq_id, self.start+1, self.end)
            else:
                return "{}:{}-{}".format(self.seq_id, self.start+1, self.end)
        else:
            if self.name:
                return "{} {}:{}-{}".format(self.name, self.seq_id, self.start, self.end)
            else:
                return "{}:{}-{}".format(self.seq_id, self.start, self.end)

    def __repr__(self):
        return "{0}:{1}-{2} ({3})".format(self.seq_id, self.start, self.end, self.strand)


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


def parse_region(location_str, base=0, name=None):
    """
    Parse a string and return a location.

    Format expected is seqid|start|end, where | is some sort of delimiter.

    If end is None, 1 based was requested.

    :param location_str: a string representing a genomic location
    :param base: a number 0 or 1
    :return: a Region
    """
    loc_match = REGEX_REGION.match(location_str.strip().replace(',', ''))
    one_base = False

    if loc_match:
        loc_groups = loc_match.groups()
        identifier = loc_groups[0]
        start_base = loc_groups[2]
        start_mult = loc_groups[3]
        end_base = loc_groups[5]
        end_mult = loc_groups[6]

        if not identifier:
            raise exceptions.G2GRegionError("Cannot parse location")

        if start_base:
            try:
                start_base = int(start_base)
            except ValueError:
                raise exceptions.G2GRegionError("Start position is not numeric")
        else:
            raise exceptions.G2GRegionError("Cannot parse start position")

        start = start_base

        if start_mult:
            if start_mult.lower() not in ['mb', 'm', 'k']:
                raise exceptions.G2GRegionError("Unknown quantifier '{0}".format(start_mult))
            start = start_base * get_multiplier(start_mult)
        else:
            if base == 1:
                start -= 1

        if end_base:
            try:
                end_base = int(end_base)
            except ValueError:
                raise exceptions.G2GRegionError("End position is not numeric")
        else:
            # special condition to allow for only 1 location to be retrieved
            # or from this start base on

            if not loc_groups[4] and not end_base and not end_mult:
                # only 1 position requested
                one_base = True
            else:
                raise exceptions.G2GRegionError("Cannot parse end position")

        if end_mult:
            if end_mult.lower() not in ['mb', 'm', 'k']:
                raise exceptions.G2GRegionError("Unknown quantifier '{0}".format(end_mult))
    else:
        raise exceptions.G2GRegionError("Cannot parse location '{0}'".format(location_str))

    if one_base:
        end = None
    else:
        end = end_base * get_multiplier(end_mult)

    LOG.debug("identifier={0}, start={1}, end={2}, base={3}, name={4}".format(identifier, start, end, base, name))

    region = Region(identifier, start, end, '+', original_base=base, name=name)

    LOG.debug("parse_region returning: region = {}".format(region))

    return region



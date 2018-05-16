# -*- coding: utf-8 -*-

#
# Collection of functions related to FASTA files
#

from __future__ import print_function

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import collections
import copy
import difflib
import os
import re
import sys
import time

import pysam

from . import exceptions
from . import gtf_db
from . import g2g
from . import g2g_utils
from . import vci

LOG = g2g.get_logger()


def filter_fai(line):
    if line.startswith("[fai"):
        return True

REGEX_FASTA_HEADER = re.compile(">(\S*)\s*(.*)")
FASTA_HEADER_FIELDS = ['id', 'description']
FastaHeader = collections.namedtuple('FastaHeader', FASTA_HEADER_FIELDS)


class FastaFile(object):
    """
    Encapsulate a Fasta file

    Delegation used with `pysam.FastaFile`
    """
    def __init__(self, filename):
        self.filename = filename
        self._fasta_file = None
        try:
            with g2g_utils.silence(filter_fai, sys.stderr):
                self._fasta_file = pysam.FastaFile(self.filename)
        except:
            self._fasta_file = pysam.FastaFile(self.filename)

        self.fai = FAI(self.filename)
        self.haploid = False
        self.diploid = False

        for record in self.fai.records:
            if record.endswith(('_L', '_R')):
                self.diploid = True
                self.haploid = False
                break

        if not self.diploid:
            self.haploid = True

    def is_diploid(self):
        return self.diploid

    def is_haploid(self):
        return self.haploid

    def fetch_list(self, reference=None, start=None, end=None, region=None):
        _start = 0 if start < 0 else start
        return list(self._fasta_file.fetch(reference=reference, start=_start, end=end, region=region))

    def __getattr__(self, name):
            return getattr(self._fasta_file, name)


class FAIEntry(object):
    def __init__(self, idx=None, length=0, offset=-1, line_length=None, line_length_bytes=None):
        """
        Initialize a `.fasta.FAIEntry` object.

        :param idx: index of fasta sequence
        :type idx: string
        :param length: length of fasta sequence
        :type length: int
        :param offset: offset in file
        :type offset: int
        :param line_length: line length
        :type line_length: int
        :param line_length_bytes: line length + newline
        :type line_length_bytes: int
        :return: `.fasta.FAIEntry`
        """
        self.idx = idx
        self.length = length
        self.offset = offset
        self.line_length = line_length
        self.line_length_bytes = line_length_bytes

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(self.idx, self.length, self.offset, self.line_length, self.line_length_bytes)


class FAI(object):
    def __init__(self, fasta_filename):
        """
        Initialize a `.fasta.FAI` file, an index into the Fasta file.

        :param file_name: the name of the Fasta file
        :type file_name: string
        :return: an initialized `.fasta.FAI` object
        """
        self.fasta_file = fasta_filename
        self.fai_file = "{0}.fai".format(fasta_filename)
        self.records = collections.OrderedDict()

        if os.path.exists(self.fai_file):
            self.read()
        else:
            self.fai_file = "{0}.gz.fai".format(fasta_filename)
            if os.path.exists(self.fai_file):
                self.read()
            else:
                raise exceptions.G2GFastaError("Unable to find fasta index file")

    def __getitem__(self, entry):
        try:
            if isinstance(entry, int):
                entry = tuple(self.records.keys())[entry]

            return self.records[entry]
        except KeyError:
            raise exceptions.G2GFastaError("{0} not in {1}.".format(entry, self.fai_file))

    def __iter__(self):
        for keys in self.records:
            yield keys

    def __repr__(self):
        return "{0}('{1}')".format(self.__class__, self.fai_file)

    def read(self):
        with open(self.fai_file) as index:
            for line in index:
                line = line.strip()
                idx, length, offset, line_length, line_length_bytes = line.split('\t')
                if idx in self.records:
                    raise ValueError("A duplicate entry was found: {0}".format(idx))
                else:
                    self.records[idx] = FAIEntry(idx, int(length), int(offset), int(line_length), int(line_length_bytes))

    def get_pos(self, seq_id, start, end):
        """
        Get the byte positions of the fasta file using the specified location

        :param seq_id: chromosome to use
        :type seq_id: string
        :param start: starting position
        :type start: int
        :param end: ending position
        :type end: int
        :return: the starting byte, ending byte, and byte length
        """
        chrom = self.records[seq_id]

        fai_entry_length = chrom.length
        fai_entry_offset = chrom.offset
        fai_entry_line_length = chrom.line_length
        fai_entry_line_length_bytes = chrom.line_length_bytes
        seq_len = end - start
        line_ratio = fai_entry_line_length * (fai_entry_line_length_bytes - fai_entry_line_length)
        newlines_total = int(fai_entry_length / line_ratio)
        newlines_before = 0

        if start > 0:
            newlines_before = int(start / line_ratio)

        newlines_to_end = int(end / line_ratio)
        byte_len_seq = newlines_to_end - newlines_before + seq_len
        byte_start = fai_entry_offset + newlines_before + start
        byte_end = fai_entry_offset + newlines_total + fai_entry_length

        return byte_start, byte_end, byte_len_seq

    def dump(self):
        print("FAI")
        print("Fasta File: {}".format(self.fasta_file))
        print("FAI File: {}".format(self.fai_file))
        for record in self.records:
            print(str(self.records[record]))


def extract(fasta_input, locations, output=None, reverse=False, complement=False, raw=False):
    """

    :param fasta_input: the name of the Fasta file
    :type fasta_input: string
    :param locations: list of `.g2g_utils.Region` objects
    :type locations: list
    :param output: name of output file, None for sys.stdout
    :type output: string
    :param reverse: reverse the extracted sequence
    :type reverse: boolean
    :param complement: complement the extracted sequence
    :type complement: boolean
    :param raw: for just the sequence to be output, no fasta line
    :type raw: boolean
    """

    start_time = time.time()

    if isinstance(fasta_input, pysam.FastaFile):
        fasta = fasta_input
    else:
        fasta_input = g2g_utils.check_file(fasta_input)
        fasta = pysam.FastaFile(fasta_input)

    LOG.debug("FASTA FILE: {0}".format(fasta.filename))

    if output:
        output = g2g_utils.check_file(output, 'w')
        fasta_out = open(output, "w")
    else:
        fasta_out = sys.stdout

    try:
        if not isinstance(locations, list):
            locations = [locations]

        for location in locations:
            LOG.debug('LOCATION: {0}'.format(location))
            if not isinstance(location, g2g.Region):
                raise exceptions.G2GRegionError("Found: {0}, which is not a Region")

            if location.seq_id not in fasta.references:
                continue

            sequence = fasta.fetch(location.seq_id, location.start, location.end)

            if location.strand == '+':
                if reverse and not complement:
                    sequence = g2g_utils.reverse_sequence(sequence)
                elif not reverse and complement:
                    sequence = g2g_utils.complement_sequence(sequence)
                elif reverse and complement:
                    sequence = g2g_utils.reverse_complement_sequence(sequence)
            else:
                sequence = g2g_utils.reverse_complement_sequence(sequence)

            if raw:
                fasta_out.write(sequence)
                fasta_out.write('\n')
            else:
                # internally storing 0 base, output 1 base
                start = location.start + 1
                end = location.end
                len_seq = len(sequence)

                # if someone requests something to long, just truncate the length
                if len_seq < end:
                    end = len_seq + start - 1

                fasta_id = ">{0}:{1}-{2}\n".format(location.seq_id, start, end)
                LOG.debug("Fasta ID: {}".format(fasta_id))
                LOG.debug("{0}...{1}".format(sequence[:10], sequence[-10:]))

                if location.name:
                    fasta_id = ">{0} {1}:{2}-{3}\n".format(location.name, location.seq_id, start, end)
                    LOG.debug("Lcation name specified, Fasta ID: {}".format(fasta_id))

                fasta_out.write(fasta_id)

                g2g_utils.write_sequence(sequence, fasta_out)

    except exceptions.G2GRegionError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except exceptions.G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    fasta_out.close()

    LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start_time, time.time())))


def extract_id(fasta_input, identifier, output=None, reverse=False, complement=False, raw=False):
    """
    """
    LOG.info("Execution start")

    start_time = time.time()

    if isinstance(fasta_input, pysam.FastaFile):
        fasta = fasta_input
    else:
        fasta_input = g2g_utils.check_file(fasta_input)
        fasta = pysam.FastaFile(fasta_input)

    LOG.debug("FASTA FILE: {0}".format(fasta.filename))

    if output:
        output = g2g_utils.check_file(output, 'w')
        fasta_out = open(output, "w")
    else:
        fasta_out = sys.stdout

    try:
        sequence = fasta[identifier]

        if reverse and not complement:
            sequence = g2g_utils.reverse_sequence(sequence)
        elif not reverse and complement:
            sequence = g2g_utils.complement_sequence(sequence)
        elif reverse and complement:
            sequence = g2g_utils.reverse_complement_sequence(sequence)

        if raw:
            fasta_out.write(sequence)
        else:
            fasta_id = ">{0}\n".format(identifier)
            fasta_out.write(fasta_id)

            for line in g2g_utils.wrap_sequence(sequence):
                fasta_out.write(line.strip())
                fasta_out.write('\n')

    except exceptions.G2GRegionError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except exceptions.G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    fasta_out.close()

    LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start_time, time.time())))


def diff_files(file1, file2):

    fasta_1 = pysam.FastaFile(file1)
    fasta_2 = pysam.FastaFile(file2)

    fa_equal = {}
    fa_diff = {}

    for seq_id1 in fasta_1.references:
        if seq_id1 in fasta_2.references:
            print("Comparing {}".format(seq_id1))

            if fasta_1.fetch(seq_id1) == fasta_2.fetch(seq_id1):
                fa_equal[seq_id1] = seq_id1
            else:
                fa_diff[seq_id1] = seq_id1

    print("# EQUAL: {0}".format(len(fa_equal)))
    print("# DIFF: {0}".format(len(fa_diff)))

    for seq in fa_diff:
        print(seq)


def diff_sequence(sequence1, sequence2):
    """

    cases=[('afrykanerskojęzyczny', 'afrykanerskojęzycznym'),
           ('afrykanerskojęzyczni', 'nieafrykanerskojęzyczni'),
           ('afrykanerskojęzycznym', 'afrykanerskojęzyczny'),
           ('nieafrykanerskojęzyczni', 'afrykanerskojęzyczni'),
           ('nieafrynerskojęzyczni', 'afrykanerskojzyczni'),
           ('abcdefg','xac')]

    for a,b in cases:
        print('{} => {}'.format(a,b))
        for i,s in enumerate(difflib.ndiff(a, b)):
            if s[0]==' ': continue
            elif s[0]=='-':
                print(u'Delete "{}" from position {}'.format(s[-1],i))
            elif s[0]=='+':
                print(u'Add "{}" to position {}'.format(s[-1],i))
        print()
    """
    for i, s in enumerate(difflib.ndiff(sequence1, sequence2)):
        print(i, s)

        if s[0]==' ':
            continue
        elif s[0] == '-':
            print(u"Delete '{}' from position {}".format(s[-1], i))
        elif s[0] == '+':
            print(u"Add '{}' to position {}".format(s[-1], i))
    print()


def compare_files(file_1, file_2, enid):
    f1 = pysam.FastaFile(file_1)
    f2 = pysam.FastaFile(file_2)

    fa_equal={}
    fa_diff={}


    c = 0

    if enid:
        pass
    else:
        for ensembl_id in f1:
            seq_1 = f1[ensembl_id]
            seq_2 = f2[ensembl_id]
            if len(str(seq_1)) == len(str(seq_2)):
            #if str(seq_1) == str(seq_2):
                fa_equal[ensembl_id] = {'seq_1':str(seq_1), 'seq_2':str(seq_2)}
            else:
                fa_diff[ensembl_id] = {'seq_1':str(seq_1), 'seq_2':str(seq_2)}
                if abs(len(str(seq_1)) - len(str(seq_2))) > 30:
                    print(ensembl_id)

    print("# EQUAL: {0}".format(len(fa_equal)))
    print("# DIFF: {0}".format(len(fa_diff)))


def get_pos(fai, chromosome, start, end):
    """
    Get the byte positions of the fasta file using the specified location

    :param fasta_index: fasta index file used for sequence retrieval
    :type fasta_index: :class:`.fasta.FAI`
    :param chromosome: chromosome to use
    :type chromosome: string
    :param start: starting position
    :type start: int
    :param end: ending position
    :type end: int
    :return: a tuple of starting byte, ending byte, and byte length
    """
    chrom = fai.records[chromosome]
    fai_entry_length = chrom.length
    fai_entry_offset = chrom.offset
    fai_entry_line_length = chrom.line_length
    fai_entry_line_length_bytes = chrom.line_length_bytes
    seq_len = end - start
    line_ratio = fai_entry_line_length * (fai_entry_line_length_bytes - fai_entry_line_length)
    newlines_total = int(fai_entry_length / line_ratio)
    newlines_before = 0
    if start > 0:
        newlines_before = int(start / line_ratio)
    newlines_to_end = int(end / line_ratio)
    byte_len_seq = newlines_to_end - newlines_before + seq_len
    byte_start = fai_entry_offset + newlines_before + start
    byte_end = fai_entry_offset + newlines_total + fai_entry_length
    return byte_start, byte_end, byte_len_seq


def parse_fasta_header(fasta_header):
    """

    :param fasta_header:
    :return:
    """
    if not fasta_header:
        return None

    match = REGEX_FASTA_HEADER.match(fasta_header.strip())

    if match:
        groups = match.groups()
        return FastaHeader(groups[0], groups[1] if groups[1] else None)

    return None



def reformat(fasta_input, output_file, length=60):
    fasta_input = g2g_utils.check_file(fasta_input)
    output = sys.stdout

    if output_file:
        output_file = g2g_utils.check_file(output_file, 'w')
        output = open(output_file, 'w')

    with open(fasta_input, 'r') as input:
        new_sequence = StringIO()
        new_header = '>UNKNOWN'

        for line in input:
            line = line.strip()

            if line[0] == '>':
                if len(new_sequence.getvalue()) > 0:
                    output.write(new_header)
                    output.write('\n')
                    g2g_utils.write_sequence(new_sequence.getvalue(), output, length)

                    LOG.info('Reformatting {}'.format(new_header))

                new_sequence = StringIO()
                new_header = line
            else:
                new_sequence.write(line)

        LOG.info('Reformatting {}'.format(new_header))
        output.write(new_header)
        output.write('\n')
        g2g_utils.write_sequence(new_sequence.getvalue(), output, length)


    output.close()


def fasta_extract_transcripts(fasta_file, database_file, output, raw=False):
    start = time.time()

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file = g2g_utils.check_file(database_file)

    fasta_out = sys.stdout

    if output:
        output = g2g_utils.check_file(output, 'w')
        fasta_out = open(output, "w")

    LOG.info("FASTA FILE: {0}".format(fasta.filename))
    LOG.info("DATABASE FILE: {0}".format(database_file))
    LOG.info("OUTPUT FILE: {0}".format(fasta_out.name))

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file)

        for i, transcript in enumerate(transcripts):
            LOG.debug("Transcript={0}".format(transcript))

            if transcript.seqid not in fasta.references:
                LOG.debug('skipping, transcript seqid not in fasta references')
                LOG.debug('{} not in {}'.format(transcript.seqid, fasta.references))
                continue

            new_sequence = StringIO()
            locations = []
            for ensembl_id, exon in transcript.exons.items():
                LOG.debug("Exon ID={0};{1}".format(ensembl_id, exon))

                partial_seq = fasta.fetch(exon.seqid, exon.start-1, exon.end)
                partial_seq_str = str(partial_seq)

                if transcript.strand == 1:
                    new_sequence.write(partial_seq_str)
                else:
                    partial_seq_str = str(g2g_utils.reverse_complement_sequence(partial_seq))
                    new_sequence.write(partial_seq_str)

                LOG.debug("{0}:{1}-{2} (Length: {3})\n{4}".format(exon.seqid, exon.start, exon.end, len(partial_seq), partial_seq_str))
                locations.append("{0}:{1}-{2}".format(exon.seqid, exon.start, exon.end))

            if raw:
                fasta_out.write(new_sequence.getvalue())
            else:
                fasta_id = ">{0} {1}|{2}\n".format(transcript.ensembl_id, '-' if transcript.strand == -1 else '+', "|".join(locations))
                fasta_out.write(fasta_id)

                for line in g2g_utils.wrap_sequence(new_sequence.getvalue()):
                    fasta_out.write(line.strip())
                    fasta_out.write('\n')

    except exceptions.G2GValueError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except exceptions.G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start, time.time())))


def fasta_extract_exons(fasta_file, database_file, output, raw=False):
    start = time.time()

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file = g2g_utils.check_file(database_file)

    fasta_out = sys.stdout

    if output:
        output = g2g_utils.check_file(output, 'w')
        fasta_out = open(output, "w")

    LOG.info("FASTA FILE: {0}".format(fasta.filename))
    LOG.info("DATABASE FILE: {0}".format(database_file))
    LOG.info("OUTPUT FILE: {0}".format(fasta_out.name))

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file)
        for i, transcript in enumerate(transcripts):

            if transcript.seqid not in fasta.references:
                continue

            for ensembl_id, exon in transcript.exons.items():
                LOG.debug("Exon={0}".format(exon))

                partial_seq = fasta.fetch(exon.seqid, exon.start-1, exon.end)
                partial_seq_str = partial_seq

                if transcript.strand == -1:
                    partial_seq_str = str(g2g_utils.reverse_complement_sequence(partial_seq))

                LOG.debug("{0}:{1}-{2} (Length: {3})\n{4}".format(exon.seqid, exon.start, exon.end, len(partial_seq), partial_seq_str))

                if raw:
                    fasta_out.write(partial_seq_str)
                else:
                    fasta_id = ">{0} {1}:{2}-{3}\n".format(exon.ensembl_id, exon.seqid, exon.start, exon.end)
                    fasta_out.write(fasta_id)

                    for line in g2g_utils.wrap_sequence(partial_seq_str):
                        fasta_out.write(line.strip())
                        fasta_out.write('\n')

    except exceptions.G2GValueError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except exceptions.G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start, time.time())))






if __name__ == '__main__':
    g2g.configure_logging(5)
    seek_and_destroy(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

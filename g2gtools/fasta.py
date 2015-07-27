# -*- coding: utf-8 -*-

#
# Collection of functions related to FASTA files
#

import difflib
import sys
import time
import os

from collections import OrderedDict


try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from .exceptions import G2GFastaError, G2GLocationError
from .g2g_utils import Location, format_time, get_logger, complement_sequence, reverse_sequence, reverse_complement_sequence, wrap_sequence
from .g2g_fileutils import check_file

import pysam

LOG = get_logger()


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
        self.index_file = "{0}.fai".format(fasta_filename)
        self.records = OrderedDict()

        if os.path.exists(self.index_file):
            self.read()
        else:
            raise G2GFastaError("Unable to find fasta index file")

    def __getitem__(self, entry):
        try:
            if isinstance(entry, int):
                entry = tuple(self.records.keys())[entry]

            return self.records[entry]
        except KeyError:
            raise G2GFastaError("{0} not in {1}.".format(entry, self.index_file))

    def __iter__(self):
        for keys in self.records:
            yield keys

    def __repr__(self):
        return "{0}('{1}')".format(self.__class__, self.index_file)

    def read(self):
        with open(self.index_file) as index:
            for line in index:
                line = line.strip()
                idx, length, offset, line_length, line_length_bytes = line.split('\t')
                if idx in self.records:
                    raise ValueError("A duplicate entry was found: {0}".format(idx))
                else:
                    self.records[idx] = FAIEntry(idx, int(length), int(offset), int(line_length), int(line_length_bytes))


def get(fasta, location):
    """
    0-based,exclusive extraction from a Fasta file.

    ACGTACGTACGTACGT

    chr1:4-9, would yield

    ----ACGTA-------
    0123456789012345

    ACGTA

    chr1:1, would yield

    -C--------------
    0123456789012345

    C

    chr1:1-2, would yield

    -C--------------
    0123456789012345

    C (Not CG, because it is exclusive).  If you wanted CG, you would request chr1:1-3.
    """

    try:
        if location.end and location.start == location.end:
            return fasta[location.seqid][location.start]
        elif not location.end:
            return fasta[location.seqid][location.start]
        else:
            return fasta[location.seqid][location.start:location.end]
    except G2GLocationError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    return None


def extract(fasta_input, locations, output=None, reverse=False, complement=False, raw=False):
    """

    :param fasta_input: the name of the Fasta file
    :type fasta_input: string
    :param locations: list of `.g2g_utils.Location` objects
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
        fasta_input = check_file(fasta_input)
        fasta = pysam.FastaFile(fasta_input)

    LOG.debug("FASTA FILE: {0}".format(fasta.filename))

    if output:
        output = check_file(output, 'w')
        fasta_out = open(output, "w")
    else:
        fasta_out = sys.stdout

    try:
        if not isinstance(locations, list):
            locations = [locations]

        for location in locations:
            LOG.debug('LOCATION: {0}'.format(location))
            if not isinstance(location, Location):
                raise G2GLocationError("Found: {0}, which is not a Location")

            if location.seqid not in fasta.references:
                continue

            sequence = fasta.fetch(location.seqid, location.start, location.end)

            if location.strand == '+':
                if reverse and not complement:
                    sequence = reverse_sequence(sequence)
                elif not reverse and complement:
                    sequence = complement_sequence(sequence)
                elif reverse and complement:
                    sequence = reverse_complement_sequence(sequence)
            else:
                sequence = reverse_complement_sequence(sequence)

            if raw:
                fasta_out.write(sequence)
                fasta_out.write('\n')
            else:
                # internally storing 0 base, outputting 1 base
                start = location.start + 1
                end = location.end

                fasta_id = ">{0} {1}-{2}\n".format(location.seqid, start, end)
                LOG.debug(fasta_id)
                LOG.debug("{0}...{1}".format(sequence[:10], sequence[-10:]))

                if location.name:
                    fasta_id = ">{0} {1}:{2}-{3}\n".format(location.name, location.seqid, start, end)

                fasta_out.write(fasta_id)

                for line in wrap_sequence(sequence):
                    fasta_out.write(line.strip())
                    fasta_out.write('\n')

    except G2GLocationError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    fasta_out.close()

    LOG.info("Execution complete: {0}".format(format_time(start_time, time.time())))


def extract_id(fasta_input, identifier, output=None, reverse=False, complement=False, raw=False):
    """
    """
    LOG.info("Execution start")

    start_time = time.time()

    if isinstance(fasta_input, pysam.FastaFile):
        fasta = fasta_input
    else:
        fasta_input = check_file(fasta_input)
        fasta = pysam.FastaFile(fasta_input)

    LOG.debug("FASTA FILE: {0}".format(fasta.filename))

    if output:
        output = check_file(output, 'w')
        fasta_out = open(output, "w")
    else:
        fasta_out = sys.stdout

    try:
        sequence = fasta[identifier]

        if reverse and not complement:
            sequence = reverse_sequence(sequence)
        elif not reverse and complement:
            sequence = complement_sequence(sequence)
        elif reverse and complement:
            sequence = reverse_complement_sequence(sequence)

        if raw:
            fasta_out.write(sequence)
        else:
            fasta_id = ">{0}\n".format(identifier)
            fasta_out.write(fasta_id)

            for line in wrap_sequence(sequence):
                fasta_out.write(line.strip())
                fasta_out.write('\n')

    except G2GLocationError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    fasta_out.close()

    LOG.info("Execution complete: {0}".format(format_time(start_time, time.time())))


def diff(sequence1, sequence2):
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
        print i, s

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
                    print ensembl_id

    print "# EQUAL: {0}".format(len(fa_equal))
    print "# DIFF: {0}".format(len(fa_diff))





if __name__ == '__main__':
    #diff(sys.argv[1], sys.argv[2])
    compare_files(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None)
    #create_compare_fasta_file(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


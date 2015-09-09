# -*- coding: utf-8 -*-

#
# Collection of functions related to chain files
#
# http://genome.ucsc.edu/goldenpath/help/chain.html
#


from collections import namedtuple
import os
import sys

from bx.intervals.intersection import Interval, IntervalTree

from .exceptions import G2GChainFileError
from .g2g_utils import configure_logging, get_logger, intersect_regions
from .g2g_fileutils import open_resource

LOG = get_logger()

# named tuple rather than class
IntervalMapping = namedtuple('IntervalMapping', ['from_chr', 'from_start', 'from_end', 'from_seq',
                                                 'to_chr', 'to_start', 'to_end', 'to_seq', 'same_bases', 'vcf_pos'])

CHAIN_STRING = "chain 1000 {from_chr} {from_length} + {from_start} {from_end} " + \
               "{to_chr} {to_length} + {to_start} {to_end} {id}"

CHAIN_HEADER_CHAIN = 0
CHAIN_HEADER_SCORE = 1
CHAIN_HEADER_TNAME = 2
CHAIN_HEADER_TSIZE = 3
CHAIN_HEADER_TSTRAND = 4
CHAIN_HEADER_TSTART = 5
CHAIN_HEADER_TEND = 6
CHAIN_HEADER_QNAME = 7
CHAIN_HEADER_QSIZE = 8
CHAIN_HEADER_QSTRAND = 9
CHAIN_HEADER_QSTART = 10
CHAIN_HEADER_QEND = 11
CHAIN_HEADER_ID = 12

class ChainEntry(object):
    def __init__(self):
        self.lines = []

    def __str__(self):
        return "\n".join(map(str, self.lines))


class ChainFile:
    """
    Encapsulate a chain file
    """
    def __init__(self, file_name, seq_ids=None, reverse=False):
        """
        Initialize a `.chain.ChainFile` object from a file.

        :param file_name: the name of the chain file
        :type file_name: string
        :param seq_ids: the ids of the sequences (chromosomes in most cases) to keep
        :type seq_ids: list
        :param reverse: `True` to parse in reverse
        :type reverse: boolean
        :return: and initialized `.chain.ChainFile` object
        """
        if not file_name:
            raise G2GChainFileError("Chain File must have a name")

        self.file_name = file_name
        self.dir, self.name = os.path.split(self.file_name)
        self.mapping_tree = None
        self.valid = False
        self.chrom_size_from = {}
        self.chrom_size_to = {}
        self.version = None
        self.headers = []
        self.seq_ids = seq_ids
        self.is_reversed = reverse

        self.parse(self.is_reversed)

    def get_seqids(self):
        if self.mapping_tree:
            return self.mapping_tree.keys()
        return None

    def parse(self, reverse):
        """
        Parse a chain file into mapping structure.

        The following documentation should explain a chain file format.::

            Header Line

            chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id

            The initial header line starts with the keyword chain, followed by 11 required attribute values,
            and ending with a blank line. The attributes include:

            score -- chain score
            tName -- chromosome (reference sequence)
            tSize -- chromosome size (reference sequence)
            tStrand -- strand (reference sequence)
            tStart -- alignment start position (reference sequence)
            tEnd -- alignment end position (reference sequence)
            qName -- chromosome (query sequence)
            qSize -- chromosome size (query sequence)
            qStrand -- strand (query sequence)
            qStart -- alignment start position (query sequence)
            qEnd -- alignment end position (query sequence)
            id -- chain ID

            The alignment start and end positions are represented as zero-based half-open intervals. For example, the
            first 100 bases of a sequence would be represented with start position = 0 and end position = 100, and the
            next 100 bases would be represented as start position = 100 and end position = 200. When the strand value
            is "-", position coordinates are listed in terms of the reverse-complemented sequence.

            Alignment Data Lines

            Alignment data lines contain three required attribute values:

            size dt dq
            size -- the size of the ungapped alignment
            dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
            dq -- the difference between the end of this block and the beginning of the next block (query sequence)

            NOTE: The last line of the alignment section contains only one number: the ungapped alignment size of
            the last block.

            Example:

            chain 4900 chrY 58368225 + 25985406 25985566 chr5 151006098 - 43549808 43549970 2
              16      0       2
              60      4       0
              10      0       4
              70

        :param file_name: the file to open and parse
        :param reverse: True to reverse
        :return: a ChainFile object
        """
        mapping_tree = {}
        chrom_size_from = {}
        chrom_size_to = {}
        chain_id_line = None
        line_no = 0
        process_seqid = False

        self.headers = []

        _chain_header_chain = CHAIN_HEADER_CHAIN
        _chain_header_score = CHAIN_HEADER_SCORE
        _chain_header_tname = CHAIN_HEADER_TNAME if not reverse else CHAIN_HEADER_QNAME
        _chain_header_tsize = CHAIN_HEADER_TSIZE if not reverse else CHAIN_HEADER_QSIZE
        _chain_header_tstrand = CHAIN_HEADER_TSTRAND if not reverse else CHAIN_HEADER_QSTRAND
        _chain_header_tstart = CHAIN_HEADER_TSTART if not reverse else CHAIN_HEADER_QSTART
        _chain_header_tend = CHAIN_HEADER_TEND if not reverse else CHAIN_HEADER_QEND
        _chain_header_qname = CHAIN_HEADER_QNAME if not reverse else CHAIN_HEADER_TNAME
        _chain_header_qsize = CHAIN_HEADER_QSIZE if not reverse else CHAIN_HEADER_TSIZE
        _chain_header_qstrand = CHAIN_HEADER_QSTRAND if not reverse else CHAIN_HEADER_TSTRAND
        _chain_header_qstart = CHAIN_HEADER_QSTART if not reverse else CHAIN_HEADER_TSTART
        _chain_header_qend = CHAIN_HEADER_QEND if not reverse else CHAIN_HEADER_TEND
        _chain_header_id = CHAIN_HEADER_ID

        _size = 0
        _dt = 1 if not reverse else 2
        _dq = 2 if not reverse else 1
        _same_bases = 3
        _dt_bases = 4 if not reverse else 5
        _dq_bases = 5 if not reverse else 4
        _vcf_pos = 6

        for line in open_resource(self.file_name):
            line = line.strip()
            line_no += 1

            if not line:
                continue

            if line.startswith(('#', ' ')):
                LOG.debug("Header: {0}".format(line))
                elems = line.split()
                if elems[0] == '#' and len(elems) > 1:
                    header_values = elems[1].strip().split('=')
                    if len(header_values) <= 1:
                        continue
                    self.headers.append((header_values[0].strip(), header_values[1].strip()))
                    if 'VERSION' == header_values[0].strip():
                        try:
                            self.version = int(header_values[1].strip())
                        except ValueError, e:
                            raise G2GChainFileError("Invalid chain file version '{0}'".format(self.version))
                continue

            fields = line.split()

            if fields[_chain_header_chain] == 'chain' and len(fields) in [12, 13]:
                # chain line
                LOG.debug(line)

                chain_id_line = line
                score = fields[_chain_header_score]

                # source
                t_name = fields[_chain_header_tname]

                if self.seq_ids:
                    if t_name in self.seq_ids:
                        LOG.debug("FOUND")
                        process_seqid = True
                    else:
                        process_seqid = False
                else:
                    process_seqid = True

                t_size = int(fields[_chain_header_tsize])
                t_strand = fields[_chain_header_tstrand]

                if t_strand != '+':
                    raise G2GChainFileError("tStrand in chain file must be '+'. Line #{0}: {1}".format(line_no, line))

                t_start = int(fields[_chain_header_tstart])
                t_end = int(fields[_chain_header_tend])

                #target
                q_name = fields[_chain_header_qname]
                q_size = int(fields[_chain_header_qsize])
                q_strand = fields[_chain_header_qstrand]

                if q_strand not in ['+', '-']:
                    raise G2GChainFileError("qStrand in chain file must be '+' or '-'. Line #{0}: {1}".format(line_no, line))

                q_start = int(fields[_chain_header_qstart])
                q_end = int(fields[_chain_header_qend])

                chrom_size_from[t_name] = t_size
                chrom_size_to[q_name] = q_size

                id = None if len(fields) == 12 else fields[12]

                if t_name not in mapping_tree:
                    mapping_tree[t_name] = IntervalTree()

                t_from, q_from = t_start, q_start

            elif fields[0] != 'chain' and (len(fields) == 3 or len(fields) == 7):
                # normal data line

                if not process_seqid:
                    continue

                if not self.version:
                    if len(fields) == 3:
                        self.version = 1
                    elif len(fields) == 7:
                        self.version = 2

                size, dt, dq = int(fields[_size]), int(fields[_dt]), int(fields[_dq])

                if len(fields) == 3:
                    if q_strand == '+':
                        mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_from, q_from+size, q_strand, None, None, None, None)))
                    elif q_strand == '-':
                        mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_size - (q_from+size), q_size - q_from, q_strand, None, None, None, None)))
                else:
                    if q_strand == '+':
                        mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_from, q_from+size, q_strand, fields[_dt_bases], fields[_dq_bases], fields[_same_bases], fields[_vcf_pos])))
                    elif q_strand == '-':
                        mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_size - (q_from+size), q_size - q_from, q_strand, fields[_dt_bases], fields[_dq_bases], fields[_same_bases], fields[_vcf_pos])))

                t_from += size + dt
                q_from += size + dq
            elif fields[0] != 'chain' and len(fields) == 1:
                # ending chain line
                if not process_seqid:
                    continue

                size = int(fields[0])

                if q_strand == '+':
                    mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_from, q_from+size, q_strand, None, None, None, None)))
                elif q_strand == '-':
                    mapping_tree[t_name].insert_interval(Interval(t_from, t_from+size, (q_name, q_size - (q_from+size), q_size - q_from, q_strand, None, None, None, None)))

                if (t_from + size) != t_end or (q_from + size) != q_end:
                    LOG.debug("{0} + {1} != {2}".format(t_from, size, t_end))
                    LOG.debug("or")
                    LOG.debug("{0} + {1} != {2}".format(q_from, size, q_end))
                    raise G2GChainFileError("Alignment blocks do not match specified block sizes. {0}".format(chain_id_line))
            else:
                raise G2GChainFileError("Unexpected line in chain file. Line #{0:,}: {1}".format(line_no, line))

        LOG.debug("PROCESSED {0:,} lines".format(line_no))
        self.valid = True
        self.mapping_tree = mapping_tree
        self.chrom_size_from = chrom_size_from
        self.chrom_size_to = chrom_size_to

    def find_mappings(self, chromosome, start, end):
        """
        Find mapping from source to target

        :param chromosome: chromosome
        :param start: start position
        :param end: end position
        :param strand: strand
        :return: list of IntervalMappings
        """
        mappings = []
        #LOG.debug("Finding mappings for {0}:{1}-{2}".format(chromosome, start, end))

        if chromosome not in self.mapping_tree:
            return None

        all_intervals = self.mapping_tree[chromosome].find(start, end)

        if len(all_intervals) == 0:
            #LOG.debug("No intervals found")
            return None
        else:
            for interval in all_intervals:
                #LOG.debug("Interval {0}={1}".format(len(mappings), interval))
                chromosome, real_start, real_end = intersect_regions(chromosome, start, end,
                                                                     chromosome, interval.start, interval.end)
                offset = abs(real_start - interval.start)
                size = abs(real_end - real_start)
                i_start = interval.value[1] + offset
                mappings.append(IntervalMapping(chromosome, real_start, real_end, interval.value[4],
                                                interval.value[0], i_start, i_start + size, interval.value[5], interval.value[6], interval.value[7]))

                #LOG.debug("Mapping {0}={1}".format(len(mappings)-1, mappings[-1]))

        return mappings


def collapse_entries(entries):
    """
    Recursive function to collapse chain entries.

    Example:

        [10, 20, 0]
        [2, 5, 0]
        [10]

        [5]

        [4, 0, 2]
        [10]

    would collapse to

        [10, 20, 0]
        [2, 5, 0]
        [19, 0, 2]
        [10]

    :param entries: array of chain file entries
    :return: the collapsed entries
    """
    if len(entries) == 1:
        LOG.debug('base case')
        return entries

    # combine first 2 entries
    entry1 = entries[0]
    entry2 = entries[1]

    last_line = entry1.lines[-1]
    first_line = entry2.lines[0]

    LOG.debug("last_line={0}".format(last_line))
    LOG.debug("first_line={0}".format(first_line))

    new_line = []
    lines = []

    if len(last_line) == 1 and len(first_line) == 1:
        new_line.append(last_line[0] + first_line[0])

        if len(entry1.lines) > 1:
            lines.extend(entry1.lines[:-1])

        lines.append(new_line)

    elif len(last_line) == 1 and len(first_line) != 1:

        new_line.append(last_line[0] + first_line[0])
        new_line.append(first_line[1])
        new_line.append(first_line[2])

        if len(entry1.lines) > 1:
            lines.extend(entry1.lines[:-1])
        lines.append(new_line)
        lines.extend(entry2.lines[1:])
    elif len(last_line) != 1 and len(first_line) == 1:
        LOG.error("ERROR")
    else:
        LOG.error("ERROR")

    new_entry = ChainEntry()
    new_entry.lines = lines

    new_entries = []
    new_entries.append(new_entry)
    new_entries.extend(entries[2:])

    return collapse_entries(new_entries)


class ChainIter(object):
    """
    Simple way to iterate chain file
    """
    def __init__(self, file_name, reverse=False):
        if not file_name:
            raise G2GChainFileError("A filename must be supplied")

        self.file_name = file_name
        self.reverse = reverse
        self.current_chain_header = None
        self.current_line = None
        self.current_record = None
        self.line_no = 0
        self.reader = open_resource(file_name)

    def reset(self):
        if self.reader:
            try:
                self.reader.close()
            except:
                pass

        self.current_chain_header = None
        self.current_line = None
        self.current_record = None
        self.line_no = 0
        self.reader = open_resource(self.file_name)


    def __iter__(self):
        return self

    def next(self):
        self.current_line = self.reader.next()
        self.line_no += 1

        while self.current_line.startswith("#"):
            self.line_no += 1
            self.current_line = self.reader.next()

        while len(self.current_line.strip()) == 0:
            self.line_no += 1
            self.current_line = self.reader.next()

        self.current_record = self.parse_line(self.current_line)

        return self.current_record

    def parse_line(self, line):
        if not line:
            raise G2GChainFileError("Empty line")

        fields = line.split()

        _chain_header_chain = CHAIN_HEADER_CHAIN
        _chain_header_score = CHAIN_HEADER_SCORE
        _chain_header_tname = CHAIN_HEADER_TNAME if not self.reverse else CHAIN_HEADER_QNAME
        _chain_header_tsize = CHAIN_HEADER_TSIZE if not self.reverse else CHAIN_HEADER_QSIZE
        _chain_header_tstrand = CHAIN_HEADER_TSTRAND if not self.reverse else CHAIN_HEADER_QSTRAND
        _chain_header_tstart = CHAIN_HEADER_TSTART if not self.reverse else CHAIN_HEADER_QSTART
        _chain_header_tend = CHAIN_HEADER_TEND if not self.reverse else CHAIN_HEADER_QEND
        _chain_header_qname = CHAIN_HEADER_QNAME if not self.reverse else CHAIN_HEADER_TNAME
        _chain_header_qsize = CHAIN_HEADER_QSIZE if not self.reverse else CHAIN_HEADER_TSIZE
        _chain_header_qstrand = CHAIN_HEADER_QSTRAND if not self.reverse else CHAIN_HEADER_TSTRAND
        _chain_header_qstart = CHAIN_HEADER_QSTART if not self.reverse else CHAIN_HEADER_TSTART
        _chain_header_qend = CHAIN_HEADER_QEND if not self.reverse else CHAIN_HEADER_TEND
        _chain_header_id = CHAIN_HEADER_ID

        _size = 0
        _dt = 1 if not self.reverse else 2
        _dq = 2 if not self.reverse else 1
        _same_bases = 3
        _dt_bases = 4 if not self.reverse else 5
        _dq_bases = 5 if not self.reverse else 4
        _vcf_pos = 6

        if fields[CHAIN_HEADER_CHAIN] == 'chain' and len(fields) in [12, 13]:
            # chain line
            LOG.debug(line)

            chain_id_line = line
            score = fields[_chain_header_score]

            # source
            t_name = fields[_chain_header_tname]
            t_size = int(fields[_chain_header_tsize])
            t_strand = fields[_chain_header_tstrand]

            if t_strand != '+':
                raise G2GChainFileError("tStrand in chain file must be '+'. Line #{0}: {1}".format(self.line_no, line))

            t_start = int(fields[_chain_header_tstart])
            t_end = int(fields[_chain_header_tend])

            #target
            q_name = fields[_chain_header_qname]
            q_size = int(fields[_chain_header_qsize])
            q_strand = fields[_chain_header_qstrand]

            if q_strand not in ['+', '-']:
                raise G2GChainFileError("qStrand in chain file must be '+' or '-'. Line #{0}: {1}".format(self.line_no, line))

            q_start = int(fields[_chain_header_qstart])
            q_end = int(fields[_chain_header_qend])

            id = None if len(fields) == 12 else fields[12]

            self.current_chain_header = [score, t_name, t_size, t_strand, t_start, t_end, q_name, q_size, q_strand, q_start, q_end, id]

            return self.current_chain_header
        elif fields[0] != 'chain' and (len(fields) == 3 or len(fields) == 7):
            # normal data line
            return [int(fields[_size]), int(fields[_dt]), int(fields[_dq]), fields[_same_bases], fields[_dt_bases], fields[_dq_bases], fields[_vcf_pos]]
        elif fields[0] != 'chain' and len(fields) == 1:
            return [int(fields[0])]



def abc(a):
    print a

if __name__ == '__main__':
    """
    #print 'a'
    #c = ChainFile(sys.argv[1], 2, '11')
    #print 'hi'
    import struct
    fd = open('akjahs', "w")
    for c in xrange(20):
        print str(c)
        for i in xrange(2000000):
            a = (str(c), i, i*2, str(c), i+1, i*2+1, "ACGGT", "AAD")
            struct.pack(a)
    fd.close()
    """
    configure_logging(5)
    for r in [False, True]:
        c = ChainFile(sys.argv[1], 2, reverse=r)
        c.mapping_tree['1'].traverse(abc)
        #print ""
        #print '1', 3000015, 3000099
        #print c.find_mappings('1', 3000015, 3000099)
        #print ""
        #print '1', 3000015, 3000099
        #print c.find_mappings('1', 3000015, 3000099)
        #print ""
        #print '1', 3000099, 3000015
        #print c.find_mappings('1', 3000099, 3000015)
        #print ""
        #print '1', 3000099, 3000015
        #print c.find_mappings('1', 3000099, 3000015)
        #print ""
        print 'done'



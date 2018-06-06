# -*- coding: utf-8 -*-

from __future__ import print_function

import collections
import multiprocessing
import os
import time

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from . import exceptions
from . import g2g
from . import g2g_utils

import pysam
from bx.intervals.intersection import Interval, IntervalTree


LOG = g2g.get_logger()

IntervalInfo = collections.namedtuple('IntervalInfo', ['chr', 'start', 'end', 'shared', 'inserted', 'deleted', 'pos'])

IntervalMapping = collections.namedtuple('IntervalMapping', ['from_chr', 'from_start', 'from_end', 'from_seq',
                                                             'to_chr', 'to_start', 'to_end', 'to_seq',
                                                             'same_bases', 'vcf_pos'])

global_mapping_tree = None

#
# VCIFile
#
# HEADERS
# ##EXECUTION_TIME
# ##INDEL_VCF
# ##INDEL_SNP
# ##STRAIN
# ##VCF_KEEP
# ##FILTER_PASSED
# ##FILTER_QUALITY
# ##DIPLOID
# ##PROCESSES
# #CHROM POS SHARE INS DEL FRAG
# DATA

def d(x):
    print(x)


def process_piece(args):
    filename = args['filename']
    contig = args['contig']
    reverse = args['reverse']

    try:
        _chrom = 0
        _pos = 1
        _shared = 2
        _deleted = 3 if not reverse else 4
        _inserted = 4 if not reverse else 3
        _fragment = 5

        num_lines_chrom = 0
        num_lines_processed = 0

        pos_from = 0
        pos_to = 0

        tree = IntervalTree()

        tabix_file = pysam.TabixFile(filename)
        iterator = tabix_file.fetch(contig, parser=pysam.asTuple())

        LOG.info("Parsing VCI, contig: {}".format(contig))

        for rec in iterator:
            num_lines_chrom += 1

            if len(rec) != 6:
                raise exceptions.G2GError("Unexpected line in VCI file: {0}".format(rec))

            if rec[2] == '.':
                continue

            """
            1 3000019 G . A 3000019
            1 3003170 A . T 3151
            1 3003197 A . G 27
            1 3003640 CG  GGGG  . 444
            1 3006790 G AA  . 3145
            1 3006834 G A . 42
            1 3007272 GC  C . 438
            1 3008489 T . ATC 1215
            """

            #LOG.debug("||".join(rec))

            fragment = int(rec[_fragment])
            deleted_bases = 0 if rec[_deleted] == '.' else len(rec[_deleted])
            inserted_bases = 0 if rec[_inserted] == '.' else len(rec[_inserted])

            #LOG.debug("pos_from={}, pos_to={}".format(pos_from, pos_to))
            #LOG.debug("Inserting interval {} - {}".format(pos_from, pos_from + fragment))
            interval = Interval(pos_from, pos_from + fragment,
                                IntervalInfo(contig, pos_to, pos_to + fragment, rec[_shared],
                                             rec[_deleted], rec[_inserted], rec[_pos]))

            #LOG.debug(interval)

            tree.insert_interval(interval)

            pos_from += (fragment + deleted_bases)
            pos_to += (fragment + inserted_bases)
            num_lines_processed += 1

        #LOG.debug("Parsed {0:,} lines for contig {1} in {2}".format(num_lines_processed, contig, g2g_utils.format_time(cotig_start_time, time.time())))
    except Exception as e:
        LOG.error(e)

    return {'tree': tree,
            'contig': contig}


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


class VCIFile:
    """
    Encapsulate a VCI (Variant Call Information) File
    """
    HEADERS = [
        'CREATION_TIME',
        'INDEL_VCF',
        'INDEL_SNP',
        'STRAIN',
        'VCF_KEEP',
        'FILTER_PASSED',
        'FILTER_QUALITY',
        'DIPLOID',
        'PROCESSES'
    ]

    def __init__(self, filename, mode='r', parser=None, index=None, encoding="ascii", seq_ids=None, reverse=False):
        if not filename:
            raise exceptions.G2GError("VCI File must have a name")

        self.filename = filename
        self.dir, self.name = os.path.split(self.filename)

        self.mapping_tree = {}
        self.headers = {}
        self.seq_ids = seq_ids
        self.is_reversed = reverse
        self.contigs = {}

        self._tabix_file = pysam.TabixFile(self.filename, mode=mode, parser=parser, index=index, encoding=encoding)

        g2g_utils.index_file(original_file=filename, file_format="vci", overwrite=False)

        self.parse_header()

    def __getattr__(self, name):
        return getattr(self._tabix_file, name)

    def fetch(self, reference=None, start=None, end=None, region=None, parser=None, multiple_iterators=False):
        # todo: potentially alter the reference value or the region value???
        return self._tabix_file.fetch(reference, start, end, region, parser, multiple_iterators)

    def is_diploid(self):
        if 'DIPLOID' in self.headers:
            return self.headers['DIPLOID'].lower() in ['true', 't', 'yes', '1']
        return False

    def is_haploid(self):
        if 'DIPLOID' in self.headers:
            if self.headers['DIPLOID'].lower() in ['true', 't', 'yes', '1']:
                return False
        return True

    def parse_header(self):
        self.headers = {}

        for header_line in self._tabix_file.header:
            header_line = g2g_utils.s(header_line)
            if header_line.startswith("##"):
                elems = header_line[2:].split('=')
                #print(elems)

                if elems[0] == 'CONTIG':
                    info = elems[1].split(':')
                    self.contigs[info[0]] = int(info[1])
                else:
                    self.headers[elems[0]] = elems[1]

    def get_seq_ids(self):
        if self.mapping_tree:
            return self.mapping_tree.keys()
        return None

    def parse(self, reverse):

        start = time.time()
        mapping_tree = {}
        try:
            total_num_lines_chrom = 0
            total_num_lines_processed = 0

            _chrom = 0
            _pos = 1
            _shared = 2
            _deleted = 3 if not reverse else 4
            _inserted = 4 if not reverse else 3
            _fragment = 5

            contigs = self._tabix_file.contigs

            if self.seq_ids:
                contigs = self.seq_ids

            # create a bx tree for the indels
            for contig in contigs:
                comtig_start_time = time.time()

                LOG.info("Parsing VCI, contig: {}".format(contig))
                num_lines_chrom = 0
                num_lines_processed = 0

                pos_from = 0
                pos_to = 0

                if contig not in mapping_tree:
                    mapping_tree[contig] = IntervalTree()

                iterator = None

                try:
                    iterator = self._tabix_file.fetch(contig, parser=pysam.asTuple())
                except:
                    pass

                if iterator is None:
                    continue

                for rec in iterator:
                    num_lines_chrom += 1
                    total_num_lines_chrom += 1

                    if len(rec) != 6:
                        raise exceptions.G2GError("Unexpected line in VCI file. Line #{0:,}: {1}".format(total_num_lines_chrom, rec))

                    if rec[2] == '.':
                        continue

                    """
                    1 3000019 G . A 3000019
                    1 3003170 A . T 3151
                    1 3003197 A . G 27
                    1 3003640 CG  GGGG  . 444
                    1 3006790 G AA  . 3145
                    1 3006834 G A . 42
                    1 3007272 GC  C . 438
                    1 3008489 T . ATC 1215
                    """

                    #LOG.debug("||".join(rec))

                    fragment = int(rec[_fragment])
                    deleted_bases = 0 if rec[_deleted] == '.' else len(rec[_deleted])
                    inserted_bases = 0 if rec[_inserted] == '.' else len(rec[_inserted])

                    #LOG.debug("pos_from={}, pos_to={}".format(pos_from, pos_to))
                    #LOG.debug("Inserting interval {} - {}".format(pos_from, pos_from + fragment))
                    interval = Interval(pos_from, pos_from + fragment,
                                        IntervalInfo(contig, pos_to, pos_to + fragment, rec[_shared],
                                                     rec[_deleted], rec[_inserted], rec[_pos]))

                    #LOG.debug(interval)

                    mapping_tree[contig].insert_interval(interval)

                    pos_from += (fragment + deleted_bases)
                    pos_to += (fragment + inserted_bases)
                    num_lines_processed += 1
                    total_num_lines_processed += 1


                if self.is_diploid():
                    interval = Interval(pos_from, self.contigs[contig[:-2]],
                                        IntervalInfo(contig, pos_to, pos_to + (self.contigs[contig[:-2]] - pos_from), None, None, None, None))
                else:
                    interval = Interval(pos_from, self.contigs[contig],
                                        IntervalInfo(contig, pos_to, pos_to + (self.contigs[contig] - pos_from), None, None, None, None))

                mapping_tree[contig].insert_interval(interval)

                LOG.debug("Parsed {0:,} lines for contig {1} in {2}".format(num_lines_processed, contig, g2g_utils.format_time(comtig_start_time, time.time())))

            self.valid = True
            self.mapping_tree = mapping_tree
        except Exception as e:
            g2g_utils._show_error()

        #LOG.info("Parsing complete: {0}".format(g2g_utils.format_time(start, time.time())))

    def add_to_tree(self, contig, tree):
        if contig not in self.mapping_tree:
            self.mapping_tree[contig] = tree

    def intersect_regions(self, chr1, start1, end1, chr2, start2, end2):
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
            raise exceptions.G2GRegionError("Start cannot be larger than end")

        return chr1, max(start1, start2), min(end1, end2)

    def find_before(self, chromosome, start):
        if chromosome not in self.mapping_tree:
            return None

        all_intervals = self.mapping_tree[chromosome].before(start, max_dist=10000000)

        if len(all_intervals) == 1:
            interval = all_intervals[0]
            #LOG.debug("Interval: {}".format(interval))
            chromosome, real_start, real_end = self.intersect_regions(chromosome, interval.start, interval.end,
                                                                      chromosome, interval.start, interval.end)
            offset = abs(real_start - interval.start)
            size = abs(real_end - real_start)
            i_start = interval.value[1] + offset
            mapping = IntervalMapping(chromosome, real_start, real_end, interval.value.inserted,
                                            interval.value.chr, i_start, i_start + size, interval.value.deleted,
                                            interval.value.shared, interval.value.pos)

            LOG.debug("Mapping: {}".format(mapping))

            return mapping

        return None


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

        if chromosome not in self.mapping_tree:
            LOG.debug("Chromosome {} not found in mapping tree".format(chromosome))
            LOG.debug("Available chromsomes are: {}".format(self.mapping_tree.keys()))
            return None
        else:
            LOG.debug("Chromosome {}, in mapping tree".format(chromosome))

        all_intervals = self.mapping_tree[chromosome].find(start, end)

        if len(all_intervals) == 0:
            LOG.debug("No intervals found")
            return None
        else:
            for interval in all_intervals:
                #LOG.debug("Interval {0}={1}".format(len(mappings), interval))
                chromosome, real_start, real_end = self.intersect_regions(chromosome, start, end,
                                                                          chromosome, interval.start, interval.end)
                offset = abs(real_start - interval.start)
                size = abs(real_end - real_start)
                i_start = interval.value[1] + offset
                mappings.append(IntervalMapping(chromosome, real_start, real_end, interval.value.inserted,
                                                interval.value.chr, i_start, i_start + size, interval.value.deleted,
                                                interval.value.shared, interval.value.pos))



                #LOG.debug("Mapping {0}={1}".format(len(mappings)-1, mappings[-1]))


        return mappings


def vci_query(vci_file, region, fasta_file):
    # ./bin/g2gtools vciquery -v data/mm/REF2CAST.vci.gz -r "1:13009000-13009800" -d
    start = time.time()

    vci_file = g2g_utils.check_file(vci_file, 'r')

    LOG.info("VCI File: {}".format(vci_file))
    LOG.info("Region: {}".format(region))

    vci_f = VCIFile(vci_file, seq_ids=[region.seq_id])
    vci_f.parse(False)

    mappings = vci_f.find_mappings(region.seq_id, region.start, region.end)

    for m in mappings:
        LOG.debug(m)

    start_pos = mappings[0].to_start
    end_pos = mappings[-1].to_end

    LOG.debug("Converted Region: {}:{}-{}".format(region.seq_id, start_pos + 1, end_pos + 1))

    for line in vci_f.fetch(reference=region.seq_id, start=start_pos, end=end_pos, parser=pysam.asTuple()):
        print(str(line))

    LOG.info("VCI parsed: {0}".format(g2g_utils.format_time(start, time.time())))

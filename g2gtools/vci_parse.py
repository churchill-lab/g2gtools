from . import exceptions
from . import g2g
from . import g2g_utils

from __future__ import print_function

import pysam
import collections
import time

import multiprocessing
from multiprocessing.managers import BaseManager

from bx.intervals.intersection import Interval, IntervalTree

IntervalInfo = collections.namedtuple('IntervalInfo', ['chr', 'start', 'end', 'shared', 'inserted', 'deleted', 'pos'])

IntervalMapping = collections.namedtuple('IntervalMapping', ['from_chr', 'from_start', 'from_end', 'from_seq',
                                                             'to_chr', 'to_start', 'to_end', 'to_seq',
                                                             'same_bases', 'vcf_pos'])

global_mapping_tree = {}
LOG = g2g.get_logger()


def d(x):
    print(x)

def _process_piece(filename_vci, contig, reverse):
    """

    :param contig:
    :param reverse:
    :return:
    """
    _chrom = 0
    _pos = 1
    _shared = 2
    _deleted = 3 if not reverse else 4
    _inserted = 4 if not reverse else 3
    _fragment = 5

    ret = {'contig':contig}#, 'tree':None}

    LOG.info("Chromosome: {}".format(contig))
    line_no = 0
    global global_mapping_tree


    tree = IntervalTree()

    try:

        if line_no % 10000 == 0:
            print("CONTIG {} {}".format(contig, line_no))

        pos_from = 0
        pos_to = 0

        tabix_file = pysam.TabixFile(filename_vci)

        iterator = None

        try:
            iterator = tabix_file.fetch(contig, parser=pysam.asTuple())
        except:
            LOG.debug("Exception for {}".format(contig))

        if iterator is None:
            LOG.debug("No iterator")
            return ret

        for rec in iterator:
            if len(rec) != 6:
                raise exceptions.G2GError("Unexpected line in G2G file. Line #{0:,}: {1}".format(line_no, rec))

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

            #LOG.debug(rec)

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

            line_no += 1
    except KeyboardInterrupt:
        raise exceptions.KeyboardInterruptError()
    except Exception as e:
        g2g_utils._show_error()
        LOG.error("OOOOOOOPS")

    #ret['tree'] = tree

    #print 'traversing tree...'
    #tree.traverse(d)

    Global = multiprocessing.Manager().Namespace()
    tt = Global.tree

    tt[contig] = tree
    Global.tree = tt


    LOG.debug("PROCESSED {0:,} lines for {1}".format(line_no, contig))

    return ret

def _wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return _process_piece(*args)

def parse_file(filename_vci, contigs=None, reverse=False):
    start = time.time()

    #map_tree = {}

    filename_vci = g2g_utils.check_file(filename_vci)
    tabix_file = pysam.TabixFile(filename_vci)


    if not contigs:
        contigs = tabix_file.contigs

    num_processes = multiprocessing.cpu_count()

    LOG.debug("VCI FILE: {0}".format(filename_vci))
    LOG.debug("CONTIGS: {0}".format(str(contigs)))
    LOG.debug("NUMBER OF PROCESSES: {0}".format(num_processes))
    LOG.debug("REVERSE: {0}".format(reverse))

    try:
        all_filename_vci = [filename_vci] * len(contigs)
        all_reverse = [reverse] * len(contigs)
        args = zip(all_filename_vci, contigs, all_reverse)
        pool = multiprocessing.Pool(num_processes)
        results = pool.map(_wrapper, args)

        # parse results
        #LOG.debug("Looping through results...")
        #for c in results:
        #    print str(c)
        #    LOG.debug(c['contig'])
        #    #map_tree[c['contig']] = c['tree']
        #    #c['tree'].traverse(d)

        LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start, time.time())))
    except KeyboardInterrupt:
        pool.terminate()
        raise exceptions.G2GError("Execution halted")
    except Exception as e:
        g2g_utils._show_error()
        raise exceptions.G2GError("Execution halted")

    global global_mapping_tree
    #map_tree = global_mapping_tree

    LOG.debug('mapping treeeeee')
    #print global_mapping_tree

    #for k, v in global_mapping_tree.items():
    #    print str(k)
    #    v.traverse(d)
    return global_mapping_tree
    #return map_tree



# -*- coding: utf-8 -*-

#
# A way to create a Chain file from a VCF file quickly
#

from collections import OrderedDict
from collections import Counter
import multiprocessing
import os
import sys
import time

import pysam
from pysam import TabixFile
from pysam import FastaFile

from .g2g_utils import format_time, get_logger
import g2g_fileutils as g2g_fu
from . import G2GValueError, G2GVCFError, G2GError, G2GChainFileError
from .vcf import parse_gt_new
from .chain import CHAIN_STRING
import traceback

LOG = get_logger()


class VCFtoChainInfo(object):
    def __init__(self):
        self.chromosome = None
        self.chromosome_length = 0
        self.number_vcf_lines = 0
        self.last_fragment_size = 0
        self.prev_next_ref_pos = 1
        self.chain_entries = []
        self.sums = [0, 0, 0]
        self.stats = {}
        self.end_length = 0
        self.output_file = None


class KeyboardInterruptError(Exception):
    pass


def _show_error():
    """
    show system errors
    """
    et, ev, tb = sys.exc_info()

    print "Error Type: %s" % et
    print "Error Value: %s" % ev
    while tb :
        co = tb.tb_frame.f_code
        filename = str(co.co_filename)
        line_no = str(traceback.tb_lineno(tb))
        print '    %s:%s' % (filename, line_no)
        tb = tb.tb_next


def update_stats(stats, reason):
    if not stats:
        stats = OrderedDict()

    if reason in stats:
        stats[reason] += 1
    else:
        stats[reason] = 1

    return stats


def process_piece(filename_vcf, chrom, chrom_length, sample_index, chain_info, diploid, passed, quality, vcf_keep, vcf_discard_file):
    ret = {'chrom': chrom, 'stats': {}, 'chain_info': chain_info}

    stats = OrderedDict()
    stats['ACCEPTED'] = 0

    if vcf_keep:
        vcf_discard = open(vcf_discard_file, "w")

    line_no = 0

    try:
        LOG.info("Processing Chromosome {0}...".format(chrom))
        tb = pysam.TabixFile(filename_vcf)

        for vcf_rec in tb.fetch(chrom, parser=pysam.asVCF()):
            line_no += 1

            try:
                gt = parse_gt_new(vcf_rec, sample_index)
            except:
                LOG.info("Unable to parse record, improper VCF file?")
                continue

            LOG.debug('\n')
            LOG.debug(vcf_rec)
            LOG.debug(gt)
            LOG.debug(vcf_rec[sample_index])

            if passed and 'PASS' not in vcf_rec.FILTER:

                LOG.debug("TOSSED: FILTERED ON PASS")
                stats = update_stats(stats, 'FILTERED ON PASS')

                if vcf_keep:
                    vcf_discard.write(vcf_rec)
                    vcf_discard.write("\n")
                continue

            elif quality and gt.fi == '0':

                # FI : Whether a sample was a Pass(1) or fail (0) based on FILTER values

                LOG.debug("TOSSED: FILTERED ON QUALITY")
                stats = update_stats(stats, 'FILTERED ON QUALITY')

                if vcf_keep:
                    vcf_discard.write(vcf_rec)
                    vcf_discard.write("\n")
                continue

            elif gt.left is None and gt.right is None:

                LOG.debug("TOSSED: NOT RELEVANT")
                stats = update_stats(stats, 'NOT RELEVANT')

                if vcf_keep:
                    vcf_discard.write(vcf_rec)
                    vcf_discard.write("\n")
                continue

            elif not diploid and gt.left != gt.right:
                # haploid or hexaploid
                # gt must be equal

                LOG.debug("TOSSED: HETEROZYGOUS")
                stats = update_stats(stats, 'HETEROZYGOUS')

                if vcf_keep:
                    vcf_discard.write(vcf_rec)
                    vcf_discard.write("\n")
                continue

            # START L AND R, ONLY R IF DIPLOID

            for ci, lr in chain_info.iteritems():
                if ci == 'left':
                    alt_seq = str(gt.left)
                else:
                    alt_seq = str(gt.right)

                if gt.ref == alt_seq:

                    LOG.debug("TOSSED, SAME AS REF")
                    lr.stats = update_stats(lr.stats, 'SAME AS REF')

                    if vcf_keep:
                        vcf_discard.write(vcf_rec)
                        vcf_discard.write("\n")
                    continue

                orig_alt_seq = alt_seq

                LOG.debug("SAMPLE: {0}".format(vcf_rec[sample_index]))
                LOG.debug("REF='{0}', ALT_L='{1}', ALT_R='{2}'. POS={3}".format(gt.ref, gt.left, gt.right, vcf_rec.pos))

                position = vcf_rec.pos + 1

                ref_seq = str(gt.ref)
                len_ref = len(ref_seq)
                len_alt = len(alt_seq)

                base_changes = len_ref - len_alt
                base_pos_diff = 0

                if position < lr.prev_next_ref_pos:
                    LOG.debug("TOSSED: CONFLICTING VCF ENTRIES: {0}".format(vcf_rec))

                    lr.stats = update_stats(lr.stats, 'CONFLICTING VCF ENTRIES')

                    if vcf_keep:
                        vcf_discard.write(vcf_rec)
                        vcf_discard.write("\n")

                    continue

                # find the position where the first base change is
                for n in xrange(min(len_ref, len_alt)):
                    if ref_seq[n] != alt_seq[n]:
                        base_pos_diff = n
                        break

                # if it is 0, take the minimum length
                if base_pos_diff == 0:
                    base_pos_diff = min(len_ref, len_alt)

                # add the base position difference
                position += base_pos_diff

                # recalculate the strings
                shared_bases = ref_seq[:base_pos_diff]
                ref_seq = ref_seq[base_pos_diff:]
                alt_seq = alt_seq[base_pos_diff:]

                dt = len(ref_seq)
                dq = len(alt_seq)

                next_ref_pos = position + len(ref_seq)
                fragment_size = position - lr.prev_next_ref_pos

                LOG.debug('           gt.ref: {0}'.format(gt.ref))
                LOG.debug('          ref_seq: {0}'.format(ref_seq))
                LOG.debug('               dt: {0}'.format(dt))
                LOG.debug('           gt.alt: {0}'.format(orig_alt_seq))
                LOG.debug('          alt_seq: {0}'.format(alt_seq))
                LOG.debug('               dq: {0}'.format(dq))
                LOG.debug('         position: {0}'.format(position))
                LOG.debug('prev_next_ref_pos: {0}'.format(lr.prev_next_ref_pos))
                LOG.debug('     next_ref_pos: {0}'.format(next_ref_pos))
                LOG.debug('    fragment_size: {0}'.format(fragment_size))
                LOG.debug('     base_changes: {0}'.format(base_changes))
                LOG.debug('    base_pos_diff: {0}'.format(base_pos_diff))
                LOG.debug('     shared_bases: {0}'.format(shared_bases))

                # fix any 0 length
                if fragment_size < 0:
                    LOG.debug("TOSSED: CONFLICTING VCF ENTRIES: {0}".format(vcf_rec))

                    lr.stats = update_stats(lr.stats, 'CONFLICTING VCF ENTRIES')

                    if vcf_keep:
                        vcf_discard.write(vcf_rec)
                        vcf_discard.write("\n")

                    continue

                if fragment_size != 0:
                    ref_str = ref_seq if ref_seq else '.'
                    alt_str = alt_seq if alt_seq else '.'
                    lr.chain_entries.append([fragment_size, len(ref_seq), len(alt_seq), shared_bases, ref_str, alt_str, vcf_rec.pos+1])
                else:
                    #
                    # THIS SHOULD NOT HAPPEN
                    #
                    raise G2GChainFileError('Unable to create chain file due to conflicting VCF entries')

                lr.stats = update_stats(lr.stats, 'ACCEPTED')

                LOG.debug(lr.chain_entries[-1])

                last_position = position
                lr.prev_next_ref_pos = next_ref_pos
                lr.sums[0] += fragment_size
                lr.sums[1] += dt
                lr.sums[2] += dq
                prev_line = vcf_rec

                #lr.prev_chrom = vcf_rec.contig
            chain_info[ci] = lr

        for ci, lr in chain_info.iteritems():
            #LOG.debug("CHROMOSOME[{0}] LENGTH = {1}".format(lr.prev_chrom, chrom_length))

            lr.chromosome = chrom
            lr.chromosome_length = chrom_length
            lr.last_fragment_size = chrom_length - lr.sums[0] - lr.sums[1]
            lr.end_length = lr.sums[0] + lr.last_fragment_size + lr.sums[2]
            lr.number_vcf_lines = line_no

            chain_info[ci] = lr

        if vcf_keep:
            vcf_discard.close()

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception, e:
        pass

    ret['stats'] = stats
    ret['chain_info'] = chain_info

    return ret


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


def log_stats(stats):
    LOG.info("STATISTICS")
    for s, stat in stats.iteritems():
        LOG.info("{0:,}\t\t{1}".format(stat, s))


def vcf2chain(input_file, fasta_file, strain, output_file, vcf_keep=False, passed=False, quality=False, diploid=False):
    """

    :param input_file:
    :param fasta_file:
    :param strain:
    :param output_file:
    :param vcf_keep:
    :param passed:
    :param quality:
    :param diploid:
    :return:
    """
    start = time.time()

    input_file = g2g_fu.check_file(input_file)
    fasta_file = g2g_fu.check_file(fasta_file)

    if not strain:
        raise G2GValueError("No strain was specified.")

    output_file = g2g_fu.check_file(output_file, 'w')
    output_file_dir = os.path.dirname(output_file)

    LOG.info("VCF FILE: {0}".format(input_file))
    LOG.info("FASTA FILE: {0}".format(fasta_file))
    LOG.info("CHAIN FILE: {0}".format(output_file))

    vcf_discard_file = None

    if vcf_keep:
        vcf_discard_file = "{0}.errors.vcf".format(os.path.basename(input_file))
        vcf_discard_file = os.path.join(output_file_dir, vcf_discard_file)
        LOG.info("VCF DISCARD FILE: {0}".format(vcf_discard_file))

    LOG.info("STRAIN: {0}".format(strain))
    LOG.info("PASS FILTER ON: {0}".format(str(passed)))
    LOG.info("QUALITY FILTER ON: {0}".format(str(quality)))
    LOG.info("DIPLOID: {0}".format(str(diploid)))

    if not isinstance(fasta_file, FastaFile):
        fasta_file = FastaFile(fasta_file)

    tb = TabixFile(input_file)
    sample_index = None

    for h in tb.header:
        if h[:6] == '#CHROM':
            try:
                elems = h.split('\t')
                samples = elems[9:]
                samples = dict(zip(samples, (x for x in xrange(len(samples)))))
                sample_index = samples[strain]
            except KeyError, ke:
                raise G2GVCFError("Unknown strain '{0}', valid strains are: {1}".format(strain, ", ".join(samples)))

    LOG.info("STRAIN SAMPLE INDEX: {0}".format(sample_index))

    LOG.info("Parsing VCF file...")

    # not all chromosomes/seqid will be processed if not in vcf file
    processed_seqids = OrderedDict()

    for seqid in tb.contigs:
        processed_seqids[seqid] = False

    left = VCFtoChainInfo()
    right = VCFtoChainInfo()

    chain_info = {}

    if diploid:
        left.output_file = g2g_fu.prepend_before_extension(output_file, 'left')
        right.output_file = g2g_fu.prepend_before_extension(output_file, 'right')
        chain_info['left'] = left
        chain_info['right'] = right

        g2g_fu.delete_file(left.output_file)
        g2g_fu.delete_file(right.output_file)
    else:
        left.output_file = output_file
        chain_info['left'] = left

        g2g_fu.delete_file(left.output_file)

    try:
        all_chrom = [c for c in fasta_file.references]
        all_chrom_length = [n for n in fasta_file.lengths]
        all_vcffiles = [input_file] * len(all_chrom)
        all_sample_index = [sample_index] * len(all_chrom)
        all_chain_info = [chain_info] * len(all_chrom)
        all_diploid = [diploid] * len(all_chrom)
        all_passed = [passed] * len(all_chrom)
        all_quality = [quality] * len(all_chrom)
        all_vcf_keep = [vcf_keep] * len(all_chrom)
        all_vcf_discard_file = [vcf_discard_file] * len(all_chrom)

        args = zip(all_vcffiles, all_chrom, all_chrom_length, all_sample_index, all_chain_info, all_diploid, all_passed,
                   all_quality, all_vcf_keep, all_vcf_discard_file)

        pool = multiprocessing.Pool(8)
        #pool = multiprocessing.Pool(1)
        results = pool.map(wrapper, args)

        # parse results
        total = 0

        # show stats
        for c in results:
            LOG.info("Chromosome: {0}".format(c['chrom']))
            for ci, lr in c['chain_info'].iteritems():
                A = Counter(c['stats'])
                B = Counter(lr.stats)
                log_stats(A+B)
            total += c['chain_info']['left'].number_vcf_lines

        # create header of chain file
        for c in results:
            for ci, lr in c['chain_info'].iteritems():
                outfile = open(lr.output_file, 'a')
                write = outfile.write
                write("# CHAINFILE\n")
                write("# CREATED={0}\n".format(time.strftime("%m-%d-%Y %H:%M:%S")))
                write("# PROGRAM=vcf2chain\n")
                write("# VCF={0}\n".format(input_file))
                write("# FASTA={0}\n".format(fasta_file.filename))

            break

        # loop through the results and dump to file
        for c in results:
            for ci, lr in c['chain_info'].iteritems():
                outfile = open(lr.output_file, 'a')
                write = outfile.write

                if lr.number_vcf_lines > 0:
                    write(CHAIN_STRING.format(CHAIN_STRING,
                                from_chr=c['chrom'], from_length=lr.chromosome_length,
                                from_start=0, from_end=lr.chromosome_length,
                                to_chr=c['chrom'], to_length=lr.end_length,
                                to_start=0, to_end=lr.sums[0] + lr.last_fragment_size + lr.sums[2], id=c['chrom']))
                    write("\n")

                    for c in lr.chain_entries:
                        write("\t".join(map(str, c)))
                        write("\n")
                    write(str(lr.last_fragment_size))
                    write("\n\n")
                    outfile.close()
                else:
                    chrom_length = fasta_file.get_reference_length(c['chrom'])
                    write(CHAIN_STRING.format(CHAIN_STRING,
                        from_chr=c['chrom'], from_length=chrom_length,
                        from_start=0, from_end=chrom_length,
                        to_chr=c['chrom'], to_length=chrom_length,
                        to_start=0, to_end=chrom_length, id=c['chrom']))

                    write("\n")
                    write(str(chrom_length))
                    write("\n\n")

        LOG.info("complete")

        if vcf_keep:
            vcf_discard_file.close()

        LOG.info("Parsed {0:,} total lines".format(total))
        LOG.info("Chain file created")
        LOG.info("Execution complete: {0}".format(format_time(start, time.time())))
    except KeyboardInterrupt:
        pool.terminate()
        raise G2GError("Execution halted")
    except Exception, e:
        _show_error()
        raise G2GError("Execution halted")


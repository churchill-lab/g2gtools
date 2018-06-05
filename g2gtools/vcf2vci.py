# -*- coding: utf-8 -*-

#
# A way to combine vcf files
#

from collections import OrderedDict
from past.builtins import xrange

import multiprocessing
import os
import time

import pysam

from . import exceptions
from . import fasta
from . import g2g
from . import g2g_utils
from . import vcf

LOG = g2g.get_logger()


class VCFFileInformation:
    def __init__(self, file_name=None, discard_file=None, sample_index=None):
        self.file_name = file_name
        self.discard_file = discard_file
        self.sample_index = sample_index
        self.lines = 0

    def __str__(self):
        return "{}, {}".format(self.file_name, self.sample_index)


class VCF2VCIInfo(object):
    def __init__(self):
        self.chromosome = None

        # list of VCFFileInformation
        self.vcf_files = []
        self.fasta_file = None

        # indel specific

        self.prev_next_ref_pos_right = 1
        self.prev_next_ref_pos_left = 1

        self.diploid = False
        self.passed = False
        self.quality = False
        self.vcf_keep = False

        self.output_file_left = None
        self.output_file_right = None

        self.stats_left = {}
        self.stats_right = {}


def walk_vcfs_together(readers, **kwargs):

    if 'vcf_record_sort_key' in kwargs:
        get_key = kwargs['vcf_record_sort_key']
    else:
        get_key = lambda r: (r.contig, r.pos)

    nexts = []

    for reader in readers:
        try:
            if reader:
               nexts.append(reader.next())
            else:
                nexts.append(None)
        except StopIteration:
            nexts.append(None)

    min_k = (None,)
    while any([r is not None for r in nexts]):
        next_idx_to_k = dict((i, get_key(r)) for i, r in enumerate(nexts) if r is not None)
        keys_with_prev_contig = [k for k in next_idx_to_k.values() if k[0] == min_k[0]]

        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)
        else:
            min_k = min(next_idx_to_k.values())

        min_k_idxs = set([i for i, k in next_idx_to_k.items() if k == min_k])
        yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

        for i in min_k_idxs:
            try:
                nexts[i] = readers[i].next()
            except StopIteration:
                nexts[i] = None


#: TODO: utilize stats
def update_stats(stats, reason):
    if not stats:
        stats = {}

    if reason in stats:
        stats[reason] += 1
    else:
        stats[reason] = 1

    return stats


def process_piece(merge_info):

    stats = {}

    try:
        output_file_left = None
        output_file_right = None

        if merge_info.output_file_left:
            output_file_left = open(merge_info.output_file_left, "w")

        if merge_info.output_file_right:
            output_file_right = open(merge_info.output_file_right, "w")

        mi = ['L']
        if merge_info.diploid:
            mi = ['L', 'R']

        LOG.info("Processing Chromosome {0}...".format(merge_info.chromosome))

        iterators = []
        discard_functions = []

        for i, file_info in enumerate(merge_info.vcf_files):
            vcf_tabix = pysam.TabixFile(file_info.file_name)
            try:
                vcf_iterator = vcf_tabix.fetch(merge_info.chromosome, parser=pysam.asVCF())
                iterators.append(vcf_iterator)
            except ValueError as ve:
                iterators.append(None)

            if file_info.discard_file:
                vcf_discard = open(file_info.discard_file, "w")

                def discard_record(rec):
                    vcf_discard.write(str(rec))
                    vcf_discard.write("\n")

                discard_functions.append(discard_record)
            else:
                discard_functions.append(lambda rec: None)

        n = 0

        line_numbers = 0
        for vcf_records in walk_vcfs_together(iterators):
            for i, vcf_record in enumerate(vcf_records):
                #LOG.debug(vcf_record)
                if vcf_record is None:
                    continue
                #LOG.debug(vcf_record.alt)
                #LOG.debug(type(vcf_record.alt))

                gt = vcf.parse_gt_tuple(vcf_record, merge_info.vcf_files[i].sample_index)
                #LOG.debug(gt)

                line_numbers = line_numbers + 1
                if gt.is_snp:
                    # snp
                    if merge_info.passed and 'PASS' not in vcf_record.filter:
                        discard_functions[i](vcf_record)

                    #LOG.debug("Processing SNP {}".format(vcf_record))
                    n += 1

                    if merge_info.quality and gt.fi == '0':
                        discard_functions[i](vcf_record)
                    elif gt.left is None or gt.right is None:
                        discard_functions[i](vcf_record)
                    else:
                        if merge_info.diploid:
                            # 0 i sthe same as REF and do not need
                            if gt.gt_left != 0:
                                output_file_left.write("{}_L\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, '.', vcf_record.ref, gt.left, '.'))
                            if gt.gt_right != 0:
                                output_file_right.write("{}_R\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, '.', vcf_record.ref, gt.right, '.'))
                        else:
                            if gt.gt_left == gt.gt_right and gt.gt_left != 0:
                                # ignore heterozygotes 0/1, 1/0, only process 0/0 and 1/1
                                #LOG.debug("ACCEPTED")

                                #LOG.debug('pos {} : ref {}, left {}, right {}'.format(vcf_snp.pos, vcf_snp.ref, gt.left, gt.right))
                                output_file_left.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, '.', vcf_record.ref, gt.left, '.'))
                else:
                    # indel
                    LOG.debug("Processing INDEL {}".format(vcf_record))

                    if merge_info.passed and 'PASS' not in vcf_record.filter:

                        LOG.debug("TOSSED: FILTERED ON PASS")
                        LOG.debug(vcf_record)
                        stats = update_stats(stats, 'FILTERED ON PASS')

                        discard_functions[i](vcf_record)
                        continue

                    elif merge_info.quality and gt.fi == '0':

                        # FI : Whether a sample was a Pass(1) or fail (0) based on FILTER values

                        LOG.debug("TOSSED: FILTERED ON QUALITY")
                        LOG.debug(vcf_record)
                        stats = update_stats(stats, 'FILTERED ON QUALITY')
                        discard_functions[i](vcf_record)
                        continue

                    elif gt.left is None and gt.right is None:

                        LOG.debug("TOSSED: NO STRAIN DATA")
                        LOG.debug(vcf_record)
                        stats = update_stats(stats, 'NO STRAIN DATA')
                        LOG.debug(i)
                        LOG.debug(type(vcf_record))
                        discard_functions[i](vcf_record)
                        continue

                    elif not merge_info.diploid and gt.left != gt.right:
                        # haploid or hexaploid
                        # gt must be equal

                        LOG.debug("TOSSED: HETEROZYGOUS")
                        LOG.debug(vcf_record)
                        stats = update_stats(stats, 'HETEROZYGOUS')
                        discard_functions[i](vcf_record)
                        continue

                    # START L AND R, ONLY R IF DIPLOID

                    for l_or_r in mi:
                        #LOG.debug("******************")
                        #LOG.debug(l_or_r)
                        lr_out = ''
                        if l_or_r == 'L':
                            #LOG.debug("->LEFT")
                            lr_out = '_L' if merge_info.diploid else ''
                            alt_seq = str(gt.left)
                            stats = merge_info.stats_left
                            output_file = output_file_left
                            prev_next_ref_pos = merge_info.prev_next_ref_pos_left
                        else:
                            #LOG.debug("->RIGHT")
                            lr_out = '_R' if merge_info.diploid else ''
                            alt_seq = str(gt.right)
                            stats = merge_info.stats_right
                            output_file = output_file_right
                            prev_next_ref_pos = merge_info.prev_next_ref_pos_right

                        LOG.debug("prev_next_ref_pos={}".format(prev_next_ref_pos))

                        if gt.ref == alt_seq:

                            LOG.debug("TOSSED, REF AND ALT ARE EQUAL")
                            LOG.debug(vcf_record)
                            stats = update_stats(stats, 'REF AND ALT ARE EQUAL')
                            discard_functions[i](vcf_record)
                            continue

                        orig_alt_seq = alt_seq

                        LOG.debug("SAMPLE: {0}".format(vcf_record[merge_info.vcf_files[i].sample_index]))
                        LOG.debug("REF='{0}', ALT_L='{1}', ALT_R='{2}'. POS={3}".format(gt.ref, gt.left, gt.right, vcf_record.pos))

                        position = vcf_record.pos + 1

                        ref_seq = str(gt.ref)
                        len_ref = len(ref_seq)
                        len_alt = len(alt_seq)

                        base_changes = len_ref - len_alt
                        base_pos_diff = 0

                        if position < prev_next_ref_pos:
                            LOG.debug("TOSSED: VCF ROLLBACK: {0}".format(vcf_record))
                            LOG.debug(vcf_record)

                            stats = update_stats(stats, 'VCF ROLLBACK')
                            discard_functions[i](vcf_record)
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
                        fragment_size = position - prev_next_ref_pos

                        '''

                        LOG.debug('           gt.ref: {0}'.format(gt.ref))
                        LOG.debug('          ref_seq: {0}'.format(ref_seq))
                        LOG.debug('               dt: {0}'.format(dt))
                        LOG.debug('           gt.alt: {0}'.format(orig_alt_seq))
                        LOG.debug('          alt_seq: {0}'.format(alt_seq))
                        LOG.debug('               dq: {0}'.format(dq))
                        LOG.debug('         position: {0}'.format(position))
                        LOG.debug('prev_next_ref_pos: {0}'.format(prev_next_ref_pos))
                        LOG.debug('     next_ref_pos: {0}'.format(next_ref_pos))
                        LOG.debug('    fragment_size: {0}'.format(fragment_size))
                        LOG.debug('     base_changes: {0}'.format(base_changes))
                        LOG.debug('    base_pos_diff: {0}'.format(base_pos_diff))
                        LOG.debug('     shared_bases: {0}'.format(shared_bases))
                        '''
                        # fix any 0 length
                        if fragment_size < 0:
                            #LOG.debug("TOSSED: FRAGMENT: {0}".format(vcf_record))

                            stats = update_stats(stats, 'FRAGMENT SIZE < 0')
                            discard_functions[i](vcf_record)
                            continue

                        if fragment_size != 0:
                            ref_str = ref_seq if ref_seq else '.'
                            alt_str = alt_seq if alt_seq else '.'
                            out = "{}{}\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, lr_out, vcf_record.pos+1, shared_bases, ref_str, alt_str, fragment_size)
                            LOG.debug(out)
                            output_file.write(out)
                        else:
                            #
                            # THIS SHOULD NOT HAPPEN
                            #
                            raise exceptions.G2GVCFError('Conflicting VCF entries')

                        stats = update_stats(stats, 'ACCEPTED')

                        if l_or_r == 'L':
                            merge_info.stats_left = stats
                            merge_info.prev_next_ref_pos_left = next_ref_pos
                            LOG.debug('setting merge_info.prev_next_ref_pos_left={}'.format(merge_info.prev_next_ref_pos_left))
                        else:
                            merge_info.stats_right = stats
                            merge_info.prev_next_ref_pos_right = next_ref_pos
                            LOG.debug('setting merge_info.prev_next_ref_pos_right={}'.format(merge_info.prev_next_ref_pos_right))


        if merge_info.output_file_left:
            output_file_left.close()

        if merge_info.output_file_right:
            output_file_right.close()

    except KeyboardInterrupt:
        raise exceptions.KeyboardInterruptError()
    except Exception as e:
        g2g_utils._show_error()
        raise Exception("Unknown exception")

    ret = {}
    ret['chrom'] = merge_info.chromosome
    ret['stats'] = stats
    ret['merge_info'] = merge_info
    ret['line_numbers'] = line_numbers

    return ret


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    LOG.debug(args)

    return process_piece(*args)


def log_stats(stats):
    LOG.info("STATISTICS")
    for s, stat in stats.items():
        LOG.info("{0:,}\t\t{1}".format(stat, s))


def create_vci_header(temp_directory, fasta_file, vcf_input_files, output_file, strain, vcf_keep, passed, quality, diploid, num_processes, bgzip):
    file = g2g_utils.gen_file_name("header", output_dir=temp_directory, extension='', append_time=False)
    with open(file, "w") as fd:
        fd.write("##CREATION_TIME={}\n".format(time.strftime("%m/%d/%Y %H:%M:%S")))

        for vcf_file in vcf_input_files:
            fd.write("##INPUT_VCF={}\n".format(vcf_file.file_name))

        fd.write("##FASTA_FILE={}\n".format(fasta_file))
        fd.write("##STRAIN={}\n".format(strain))
        fd.write("##VCF_KEEP={}\n".format(vcf_keep))
        fd.write("##FILTER_PASSED={}\n".format(passed))
        fd.write("##FILTER_QUALITY={}\n".format(quality))
        fd.write("##DIPLOID={}\n".format(diploid))
        fd.write("##PROCESSES={}\n".format(num_processes))

        fasta_file = fasta.FastaFile(fasta_file)

        for c in fasta_file.references:
            fd.write("##CONTIG={}:{}\n".format(c, fasta_file.get_reference_length(c)))

        fd.write("#CHROM\tPOS\tANCHOR\tINS\tDEL\tFRAG\n")
    return file


def process(vcf_files, fasta_file, output_file, strain, vcf_keep=False, passed=False, quality=False, diploid=False, num_processes=None, bgzip=False):
    start = time.time()

    output_file = g2g_utils.check_file(output_file, 'w')
    output_file_dir = os.path.dirname(output_file)

    vcf_file_inputs = []

    if vcf_files:
        for file_name in vcf_files:
            vcf_file = g2g_utils.check_file(file_name)
            LOG.info("VCF file: {0}".format(vcf_file))
            LOG.info("Checking for index file, creating if needed...")
            g2g_utils.index_file(original_file=vcf_file, file_format="vcf", overwrite=False)

            vcf_discard_file = None
            if vcf_keep:
                vcf_discard_file = "{0}.errors.vcf".format(os.path.basename(output_file))
                vcf_discard_file = os.path.join(output_file_dir, vcf_discard_file)
                LOG.info("VCF indel discard file: {0}".format(vcf_discard_file))

            vcf_file_inputs.append(VCFFileInformation(vcf_file, vcf_discard_file))

    if len(vcf_file_inputs) == 0:
        raise exceptions.G2GValueError("No VCF files.")

    if not fasta_file:
        raise exceptions.G2GValueError("No fasta file was specified.")

    if not strain:
        raise exceptions.G2GValueError("No strain was specified.")

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    LOG.info("Fasta File: {0}".format(output_file))
    LOG.info("Strain: {0}".format(strain))
    LOG.info("Pass filter on: {0}".format(str(passed)))
    LOG.info("Quality filter on: {0}".format(str(quality)))
    LOG.info("Diploid: {0}".format(str(diploid)))
    LOG.info("Number of processes: {0}".format(num_processes))
    LOG.info("Output VCI File: {0}".format(output_file))

    # not all chromosomes/seqid will be processed if not in vcf file
    processed_seq_ids = {}

    temp_directory = g2g_utils.create_temp_dir('vcf2vci', dir='.')
    LOG.debug("Temp directory: {}".format(temp_directory))

    header_file = create_vci_header(temp_directory, fasta_file, vcf_file_inputs, output_file, strain, vcf_keep, passed, quality, diploid, num_processes, bgzip)

    for i, vcf_file in enumerate(vcf_file_inputs):
        tb_file = pysam.TabixFile(vcf_file.file_name)
        for h in tb_file.header:
            h = g2g_utils.s(h)
            if h[:6] == '#CHROM':
                try:
                    elems = h.split('\t')
                    samples = elems[9:]
                    samples = dict(zip(samples, (x for x in xrange(len(samples)))))
                    vcf_file_inputs[i].sample_index = samples[strain]
                except KeyError as ke:
                    raise exceptions.G2GVCFError("Unknown strain '{0}', valid strains are: {1}".format(strain, ", ".join(samples)))

        for seq_id in tb_file.contigs:
            processed_seq_ids[seq_id] = False

    tmp_processed_seq_ids = OrderedDict()
    for k in g2g_utils.natsorted(processed_seq_ids.keys()):
        tmp_processed_seq_ids[k] = False
    processed_seq_ids = tmp_processed_seq_ids

    all_merge_info = []

    try:
        for c in processed_seq_ids:
            merge_info = VCF2VCIInfo()
            merge_info.chromosome = c

            merge_info.vcf_files = vcf_file_inputs
            merge_info.fasta_file = fasta_file

            merge_info.diploid = diploid
            merge_info.passed = passed
            merge_info.quality = quality
            merge_info.vcf_keep = vcf_keep

            if diploid:
                merge_info.output_file_left = g2g_utils.gen_file_name("chr{}.left".format(c), output_dir=temp_directory, extension='vci', append_time=False)
                merge_info.output_file_right = g2g_utils.gen_file_name("chr{}.right".format(c), output_dir=temp_directory, extension='vci', append_time=False)

                g2g_utils.delete_file(merge_info.output_file_left)
                g2g_utils.delete_file(merge_info.output_file_right)
            else:
                merge_info.output_file_left = g2g_utils.gen_file_name("chr{}.right".format(c), output_dir=temp_directory, extension='vci', append_time=False)

                g2g_utils.delete_file(merge_info.output_file_left)

            all_merge_info.append(merge_info)

        LOG.info("Parsing VCF files...")

        args = zip(all_merge_info)
        pool = multiprocessing.Pool(num_processes)
        results = pool.map(wrapper, args)

        # parse results
        total = 0
        for r in results:
            total += r['line_numbers']

        # show stats

        #TODO: make sure stats are good and show statistics

        LOG.debug("Combining temp files...")
        LOG.info("Finalizing VCI file...")

        files = [header_file]
        if diploid:
            for mi in all_merge_info:
                files.append(mi.output_file_left)
                files.append(mi.output_file_right)

            g2g_utils.concatenate_files(files, output_file, True)

            if bgzip:
                g2g_utils.bgzip_and_index_file(output_file, output_file + ".gz", delete_original=True, file_format="vcf")
        else:
            for mi in all_merge_info:
                files.append(mi.output_file_left)

            g2g_utils.concatenate_files(files, output_file, True)

            if bgzip:
                g2g_utils.bgzip_and_index_file(output_file, output_file + ".gz", delete_original=True, file_format="vcf")

#        if vcf_keep:
#            vcf_discard_file.close()

        # TODO: make sure stats are good and show statistics
        LOG.info("Parsed {0:,} total lines".format(total))
    except exceptions.KeyboardInterruptError:
        pool.terminate()
        raise exceptions.G2GError("Keyboard quit consumed")
    except KeyboardInterrupt:
        pool.terminate()
        raise exceptions.G2GError("Execution halted")
    except Exception as e:
        g2g_utils._show_error()
        raise exceptions.G2GError("Execution halted unknown error")
    finally:
        g2g_utils.delete_dir(temp_directory)
        LOG.info("VCI creation complete: {0}".format(g2g_utils.format_time(start, time.time())))


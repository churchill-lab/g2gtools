
# -*- coding: utf-8 -*-

#
# A way to patch a Fasta file quickly
#

import mmap
import multiprocessing
import os
import shutil
import time

import pysam

from .g2g_utils import format_time, get_logger
import g2g_fileutils as g2g_fu
from . import G2GValueError, G2GVCFError, G2GError
from .vcf import VCF, parse_gt, parse_vcf_line
from .fasta import FAI

LOG = get_logger()

NUM_CHUNKS = 10

# globals, really... yup :-(
fasta_index = None
mm_l = None
mm_r = None




class KeyboardInterruptError(Exception): pass






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


def process_piece(filename_vcf, chrom, start, end, sample_index, diploid, pass_only, quality):
    """
    Process a section of the VCF file specified by the location.

    :param filename_vcf: name of the VCF file
    :type filename_vcf: name of the VCF file
    :param chrom: the chromosome
    :type chrom: string
    :param start: starting position
    :type start: int
    :param end: ending position
    :type end: int
    :param sample_index: which sample index to process
    :type sample_index: int
    :param diploid: `True` to process diploid
    :type diploid: boolean
    :param pass_only: `True` to only process 'PASS' VCF records
    :type pass_only: boolean
    :param quality:
    :type quality: boolean
    :return: a dict of 2 elements, 'chrom' and 'count' specifying how many SNPs were patched
    """
    ret = {'count': 0, 'chrom': chrom}

    try:
        query = "{0}:{1}-{2}".format(chrom, start, end)
        LOG.info("Processing {0}...".format(query))
        tb = pysam.TabixFile(filename_vcf)

        global fasta_index
        for d in tb.fetch(chrom, start, end):
            rec = parse_vcf_line(d)

            if pass_only and 'PASS' not in rec.FILTER:
                continue

            byte_start, byte_end, byte_len_seq = get_pos(fasta_index, rec.CHROM, rec.POS-1, rec.POS)
            gt = parse_gt(rec, sample_index)

            #LOG.debug(rec)
            #LOG.debug(gt)

            # FI : Whether a sample was a Pass(1) or fail (0) based on FILTER values
            if quality and gt.fi == '0':
                continue

            global mm_l
            global mm_r
            if gt.left and gt.right:
                if mm_l[byte_start] != gt.ref:
                    LOG.warn("Reference at {0}:{1} has base '{2}', but VCF indicates '{3}', skipping".format(rec.CHROM, rec.POS, mm_l[byte_start], gt.ref))
                    continue

                if diploid:
                    mm_l[byte_start] = gt.left
                    mm_r[byte_start] = gt.right
                    ret['count'] += 1
                    # LOG.debug("Patching {0}:{1} ({2}) from {3} to {4},{5}".format(chrom, rec.POS-1, from_base, gt.ref, gt.left, gt.right))
                else:
                    if gt.left == gt.right:
                        # ignore heterozygotes 0/1, 1/0, only process 0/0 and 1/1
                        mm_l[byte_start] = gt.left
                        ret['count'] += 1
                        #LOG.debug("Patching {0}:{1} from {2} to {3}:{4}".format(chrom, rec.POS-1, gt.ref, gt.left, gt))


    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except ValueError, te:
        LOG.debug("No SNPS found in region {0}".format(query))
    except Exception, e:
        LOG.error('-'*80)
        LOG.error(e)
        LOG.debug("{0}, {1}, {2}".format(byte_start, byte_end, byte_len_seq))
        LOG.debug(rec)
        LOG.debug(gt)

    return ret


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


def patch(filename_original_fasta, filename_vcf, strain, filename_new_fasta_l, filename_new_fasta_r,
          num_processes, pass_only, quality, diploid):
    """
    Perform patch via multiprocess.

    :param filename_vcf: name of the VCF file
    :type filename_vcf: string
    :param strain: name of strain to use in VCF file
    :type strain: string
    :param filename_original_fasta: name of the input Fasta file
    :type filename_original_fasta: string
    :param filename_new_fasta_l: name of the output Fasta file
    :type filename_new_fasta_l: string
    :param filename_new_fasta_r: name of the output Fasta file
    :type filename_new_fasta_r: string
    :param num_processes: the number of processes to spawn
    :type num_processes: int
    :param pass_only: Only process those VCF records with a 'PASS'
    :type pass_only: boolean
    :param quality: filter on quality, FI=PASS
    :type quality: boolean
    :param diploid: don't ignore hets and create 2 files
    :type diploid: boolean
    :return: Nothing
    """

    vcf_file = VCF(filename_vcf)

    # the temporary fasta index
    global fasta_index
    fasta_index = FAI(filename_new_fasta_l)

    if strain not in vcf_file.samples:
        raise G2GVCFError("Unknown strain '{0}', valid strains are: {1}".format(strain, vcf_file.samples.keys()))

    sample_index = vcf_file.get_sample_index(strain)

    fd_l = open(filename_new_fasta_l, "r+b")
    global mm_l
    mm_l = mmap.mmap(fd_l.fileno(), 0)

    if diploid:
        fd_r = open(filename_new_fasta_r, "r+b")
        global mm_r
        mm_r = mmap.mmap(fd_r.fileno(), 0)

    all_start_pos = []
    all_end_pos = []
    all_chrom = []

    try:
        for f in fasta_index.records:
            chrom = f
            chrom_length = fasta_index.records[f].length

            # we don't want to handle to much data at a time, so chunk
            initial_chunks = range(1, chrom_length, chrom_length/NUM_CHUNKS)

            start_pos = sorted(set([i for i in initial_chunks]))
            end_pos = [i-1 for i in start_pos][1:] + [chrom_length]
            chroms = [chrom] * len(start_pos)

            all_start_pos.extend(start_pos)
            all_end_pos.extend(end_pos)
            all_chrom.extend(chroms)

        all_sample_index = [sample_index] * len(all_start_pos)
        all_diploids = [diploid] * len(all_start_pos)
        all_vcffiles = [filename_vcf] * len(all_start_pos)
        all_passonly = [pass_only] * len(all_start_pos)
        all_quality = [quality] * len(all_start_pos)

        args = zip(all_vcffiles, all_chrom, all_start_pos, all_end_pos,
                   all_sample_index, all_diploids, all_passonly, all_quality)

        pool = multiprocessing.Pool(num_processes)
        results = pool.map(wrapper, args)

        # parse results
        total = 0
        chroms = {}

        mm_l.flush()
        mm_l.close()

        if diploid:
            mm_r.flush()
            mm_r.close()

        for c in results:
            current = chroms.get(c['chrom'], 0)
            chroms[c['chrom']] = current + c['count']
            total += c['count']

        for f in fasta_index.records:
            LOG.info("Patched {0:,} SNPs on chromosome {1}".format(chroms.get(f, 0), f))

        LOG.info("Patched {0:,} SNPs total".format(total))
    except KeyboardInterrupt:
        pool.terminate()
        raise G2GError("Execution halted")
    except Exception, e:
        raise G2GError("Execution halted")


def prepare_fasta_patch(filename_fasta, filename_output, bgzip=False, diploid=False):
    """
    Initialize fasta_patch variables

    :param filename_fasta:
    :param filename_vcf:
    :param strain:
    :param filename_output:
    :param bgzip:
    :param diploid:
    :return:
    """

    filename_output = g2g_fu.check_file(filename_output, 'w')
    output_file_dir = os.path.abspath(os.path.dirname(filename_output))

    new_filename_output = filename_output

    # let's figure out what our output names will be
    if filename_output.lower().endswith('.gz'):
        # strip off .gz
        new_filename_output = filename_output[:-3]

    if not filename_output.lower().endswith('.fa'):
        raise G2GValueError("Expecting output filename extension to be either '.fa.gz' or '.fa'")


    if diploid:
        filename_output_l = g2g_fu.prepend_before_extension(new_filename_output, 'l')
        filename_output_r = g2g_fu.prepend_before_extension(new_filename_output, 'r')

        g2g_fu.delete_index_files(filename_output_l)
        g2g_fu.delete_index_files(filename_output_r)
    else:
        filename_output_l = new_filename_output
        filename_output_r = None

        g2g_fu.delete_index_files(filename_output_l)

    # at this point we are hoping for a .fa extension

    # let's figure out our input and process accordingly
    if filename_fasta.lower().endswith('.fa.gz'):
        # decompress the fasta file if it is compressed

        LOG.info("Copying and decompressing fasta file")

        # copy file and preserve gz extension for bgzip -d to work
        tmp_file_name = os.path.basename(filename_fasta)                        # something.gz
        LOG.debug("tmp_file_name={0}".format(tmp_file_name))

        tmp_fasta = os.path.join(output_file_dir, tmp_file_name)                # /path/something.fa.gz
        LOG.debug("tmp_fasta={0}".format(tmp_fasta))

        LOG.debug("COPYING {0} to {1}".format(filename_fasta, tmp_fasta))
        shutil.copy(filename_fasta, tmp_fasta)  # cp /original/something.fa.gz /output/something.fa.gz

        LOG.debug("DECOMPRESSING {0}".format(tmp_fasta))
        g2g_fu.bgzip_decompress(tmp_fasta)

        tmp_fasta = tmp_fasta[:-3]         # /path/something.fa
        LOG.debug("tmp_fasta={0}".format(tmp_fasta))

        LOG.debug("Moving '{0}' to '{1}'...".format(tmp_fasta, filename_output_l))
        shutil.move(tmp_fasta, filename_output_l)

    elif filename_fasta.lower().endswith('.fa'):
        LOG.debug("File is not compressed")

        LOG.debug("COPYING {0} to {1}".format(filename_fasta, filename_output_l))
        shutil.copy(filename_fasta, filename_output_l)
    else:
        raise G2GValueError("Expecting input filename extension to be either '.fa.gz' or '.fa'")

    if diploid:
        LOG.debug("Copying '{0}' to '{1}'...".format(filename_output_l, filename_output_r))
        shutil.copy(filename_output_l, filename_output_r)

    # build a temporary fasta index
    pysam.FastaFile(filename_output_l)

    return filename_output_l, filename_output_r



def fasta_patch(filename_fasta, filename_vcf, strain, filename_output, bgzip=False,
                num_processes=None, pass_only=False, quality=False, diploid=False):
    """
    Patch a Fasta file by replacing the bases where the SNPs are located in the VCF file.

    :param filename_fasta: name of the input Fasta file
    :type filename_fasta: string
    :param filename_vcf: name of the VCF file
    :type filename_vcf: string
    :param strain: name of strain to use in VCF file
    :type strain: string
    :param filename_output: name of the output Fasta file
    :type filename_output: string
    :param bgzip: compress file in BGZIP format
    :type bgzip: boolean
    :param num_processes: the number of processes to spawn
    :type num_processes: int
    :param pass_only: Only process those VCF records with a 'PASS'
    :type pass_only: boolean
    :param quality: filter on quality, FI=PASS
    :type quality: boolean
    :param diploid: don't ignore hets and create 2 files
    :type diploid: boolean
    :return: Nothing
    """
    start = time.time()

    filename_fasta = g2g_fu.check_file(filename_fasta)
    filename_vcf = g2g_fu.check_file(filename_vcf)

    LOG.info("INPUT FASTA FILE: {0}".format(filename_fasta))
    LOG.info("VCF FILE: {0}".format(filename_vcf))
    LOG.info("STRAIN: {0}".format(strain))
    LOG.info("PASS FILTER ON: {0}".format(str(pass_only)))
    LOG.info("QUALITY FILTER ON: {0}".format(str(quality)))
    LOG.info("DIPLOID: {0}".format(str(diploid)))

    if not strain:
        raise G2GValueError("No strain was specified.")

    filename_output_l, filename_output_r = prepare_fasta_patch(filename_fasta, filename_output, bgzip, diploid)

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    LOG.info("NUMBER OF PROCESSES: {0}".format(num_processes))
    if bgzip:
        if diploid:
            LOG.info("OUTPUT FASTA FILES: {0}.gz".format(filename_output_l))
            LOG.info("                    {0}.gz".format(filename_output_r))
        else:
            LOG.info("OUTPUT FASTA FILE: {0}.gz".format(filename_output_l))
    else:
        if diploid:
            LOG.info("OUTPUT FASTA FILES: {0}".format(filename_output_l))
            LOG.info("                    {0}".format(filename_output_r))
        else:
            LOG.info("OUTPUT FASTA FILE: {0}".format(filename_output_l))

    LOG.info("Patching...")

    try:
        patch(filename_fasta, filename_vcf, strain, filename_output_l, filename_output_r,
              num_processes, pass_only, quality, diploid)

        LOG.info("Patching complete")

        # remove the fai
        LOG.debug("removing the FAI index for {0}".format(g2g_fu.delete_index_files(filename_output_l)))
        g2g_fu.delete_index_files(filename_output_l)

        # move temp to final destination
        if bgzip:
            LOG.info("Compressing and indexing...")
            g2g_fu.bgzip_index(filename_output_l, "{0}.gz".format(filename_output_l), 'fa')
            if diploid:
                g2g_fu.bgzip_index(filename_output_r, "{0}.gz".format(filename_output_r), 'fa')

        LOG.info("Execution complete: {0}".format(format_time(start, time.time())))
    except Exception, e:
        LOG.debug(e)
        raise G2GError("")




# -*- coding: utf-8 -*-

#
# A way to transform a Fasta file quickly
#

import copy
import mmap
import multiprocessing
import os
import sys
import time

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import pysam

from . import compat
from . import exceptions
from . import fasta
from . import g2g
from . import g2g_utils
from . import vci

LOG = g2g.get_logger()


class FastaPatchParams(object):
    def __init__(self):
        # input region, must match the format of the input fasta file
        # if fasta file is haploid, format is chr:start-end
        # if fasta file is diploid, format is chr(_L|_R):start-end
        self.input_region = None

        # original input fasta file
        self.input_file = None

        # temporary directory
        self.temp_dir = None

        self.output_file = None
        self.output_region = None
        self.output_header = None

        # vci information
        self.vci_file = None
        self.vci_query = None
        self.reverse = False

        self.gen_temp = True

        # offset
        self.offset = 0

    def __str__(self):
        return "Input: {}\nOutput: {}\nLocation: {}\nOffset: {}".format(self.input_fasta_file, self.output_fasta_file, self.input_region, self.offset)


class FastaPatchResult(object):
    def __init__(self):
        self.output_file = None
        self.count = 0

    def __str__(self):
        return "File: {}".format(self.output_file)


def process_piece(fasta_patch_params):
    """

    :param fasta_patch_params:
    :return:
    """

    LOG.debug('fasta_patch_params.input_region = {}'.format(fasta_patch_params.input_region))
    LOG.debug('fasta_patch_params.input_file = {}'.format(fasta_patch_params.input_file))
    LOG.debug('fasta_patch_params.temp_dir = {}'.format(fasta_patch_params.temp_dir))
    LOG.debug('fasta_patch_params.output_file = {}'.format(fasta_patch_params.output_file))
    LOG.debug('fasta_patch_params.output_region = {}'.format(fasta_patch_params.output_region))
    LOG.debug('fasta_patch_params.output_header = {}'.format(fasta_patch_params.output_header))
    LOG.debug('fasta_patch_params.vci_file = {}'.format(fasta_patch_params.vci_file))
    LOG.debug('fasta_patch_params.vci_query = {}'.format(fasta_patch_params.vci_query))
    LOG.debug('fasta_patch_params.reverse = {}'.format(fasta_patch_params.reverse))
    LOG.debug('fasta_patch_params.gen_temp = {}'.format(fasta_patch_params.gen_temp))
    LOG.debug('fasta_patch_params.offset = {}'.format(fasta_patch_params.offset))

    try:
        fasta_patch_result = FastaPatchResult()
        fasta_patch_result.output_file = fasta_patch_params.output_file

        fasta_file = fasta.FastaFile(fasta_patch_params.input_file)
        vci_file = vci.VCIFile(fasta_patch_params.vci_file)

        if fasta_patch_params.gen_temp:
            tmp_fasta = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(fasta_patch_params.input_region)+'_', extension="fa", append_time=False, output_dir=fasta_patch_params.temp_dir)
            LOG.debug('tmp_fasta={}'.format(tmp_fasta))

            working = open(tmp_fasta, "w")
            if fasta_patch_params.output_header.description:
                fasta_header_string = ">{} {}".format(fasta_patch_params.output_header.id, fasta_patch_params.output_header.description)
            else:
                fasta_header_string = ">{}".format(fasta_patch_params.output_header.id)
            working.write("{}\n".format(fasta_header_string))

            LOG.debug("Fasta Fetch {}".format(fasta_patch_params.input_region))
            LOG.debug("Fasta Fetch {}".format(fasta_patch_params.input_region.start))
            sequence = fasta_file.fetch(fasta_patch_params.input_region.seq_id, fasta_patch_params.input_region.start - 1, fasta_patch_params.input_region.end)
            if len(sequence) < 50:
                LOG.debug("Fasta Fetch = {}".format(sequence))
            else:
                LOG.debug("Fasta Fetch = {} ... {}".format(sequence[:25], sequence[-25:]))
            g2g_utils.write_sequence(sequence, working)
            working.close()

            fasta.FastaFile(tmp_fasta)
            fasta_patch_result.output_file = tmp_fasta
            fasta_index = fasta.FAI(tmp_fasta)
        else:
            fasta_index = fasta.FAI(fasta_patch_params.input_file)
            tmp_fasta = fasta_patch_params.input_file

        offset = fasta_patch_params.offset
        reverse = fasta_patch_params.reverse

        byte_start, byte_end, byte_len_seq = 0,0,0
        LOG.info("Processing {0}...".format(fasta_patch_params.output_header.id))

        fd_l = open(tmp_fasta, "r+b")
        mm = mmap.mmap(fd_l.fileno(), 0)

        LOG.debug("VCI Fetch {}".format(fasta_patch_params.vci_query))
        for line in vci_file.fetch(fasta_patch_params.vci_query, parser=pysam.asTuple()):

            if line[5] != '.':
                continue

            fasta_patch_result.count += 1

            position = int(line[1]) - offset
            deleted_bases = line[3 if not reverse else 4]
            inserted_bases = line[4 if not reverse else 3]

            LOG.debug('int(line[1])={}'.format(int(line[1])))
            LOG.debug('position={}'.format(position))
            LOG.debug('offset={}'.format(offset))

            byte_start, byte_end, byte_len_seq = fasta_index.get_pos(fasta_patch_params.output_region.seq_id, position - 1, position)

            if compat.is_py2:
                LOG.debug("LINE: {}".format(line))
                LOG.debug('WAS {}'.format(mm[byte_start]))
                LOG.debug("Patching {0}:{1} from {2} to {3}".format(fasta_patch_params.output_region.seq_id, position - 1, deleted_bases, inserted_bases))
                mm[byte_start] = inserted_bases
            else:
                LOG.debug("LINE: {}".format(line))
                LOG.debug('WAS {}'.format(chr(mm[byte_start])))
                LOG.debug("Patching {0}:{1} from {2} to {3}".format(fasta_patch_params.output_region.seq_id, position - 1, deleted_bases, inserted_bases))
                mm[byte_start] = ord(inserted_bases)


    except KeyboardInterrupt:
        raise exceptions.KeyboardInterruptError()
    except ValueError as te:
        LOG.debug("No SNPS found in region {0}".format(fasta_patch_params.vci_query))
    except Exception as e:
        g2g_utils._show_error()
        LOG.debug("{0}, {1}, {2}".format(byte_start, byte_end, byte_len_seq))
        LOG.debug(line)

    mm.flush()
    mm.close()
    return fasta_patch_result

    #LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start_time, time.time())))


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


def prepare_fasta_patch(filename_fasta, filename_output):
    """
    Initialize fasta_patch variables

    """
    filename_output = g2g_utils.check_file(filename_output, 'w')
    output_file_dir = os.path.abspath(os.path.dirname(filename_output))

    new_filename_output = filename_output

    # let's figure out what our output name will be
    if new_filename_output.lower().endswith('.gz'):
        # strip off .gz
        new_filename_output = new_filename_output[:-3]

    if not filename_output.lower().endswith('.fa'):
        raise exceptions.G2GValueError("Expecting output filename extension to be either '.fa.gz' or '.fa'")

    g2g_utils.delete_index_files(new_filename_output)

    return new_filename_output


def process(filename_fasta, filename_vci, regions, filename_output=None, bgzip=False, reverse=False, num_processes=None):
    """
    Patch a Fasta file by replacing the bases where the SNPs are located in the VCF file.

    :param filename_fasta: name of the input Fasta file
    :type filename_fasta: string
    :param filename_g2g: name of the G2G file
    :type filename_g2g: string
    :param filename_output: name of the output Fasta file
    :type filename_output: string
    :param bgzip: compress file in BGZIP format
    :type bgzip: boolean
    :param reverse: reverse the G2G file
    :type reverse: boolean
    :param num_processes: the number of processes to spawn
    :type num_processes: int
    :return: Nothing
    """
    start = time.time()

    filename_fasta = g2g_utils.check_file(filename_fasta)
    filename_vci = g2g_utils.check_file(filename_vci)

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    LOG.info("Input Fasta File: {0}".format(filename_fasta))
    LOG.info("Input VCI File: {0}".format(filename_vci))
    LOG.info("Processes: {0}".format(num_processes))

    dump_fasta = False

    temp_directory = g2g_utils.create_temp_dir('patch_', dir='.')
    LOG.debug("Temp directory: {}".format(temp_directory))

    try:
        if filename_output:
            filename_output = g2g_utils.check_file(filename_output, 'w')

            if not regions:
                filename_output = prepare_fasta_patch(filename_fasta, filename_output)
                LOG.info("Output Fasta File: {0}".format(filename_output))
            else:
                if bgzip:
                    if filename_output.lower().endswith((".fa", ".fasta")):
                        LOG.info("Output Fasta File: {0}.gz".format(filename_output))
                    elif filename_output.lower().endswith(".gz"):
                        LOG.info("Output Fasta File: {0}".format(filename_output))
                        filename_output = filename_output[:-3]
                else:
                    LOG.info("Output Fasta File: {0}".format(filename_output))
        else:
            filename_output = g2g_utils.gen_file_name(extension="fa", append_time=False, output_dir=temp_directory)
            dump_fasta = True
            LOG.debug("Temporary fasta file: {}".format(filename_output))

        fasta_file = fasta.FastaFile(filename_fasta)
        vci_file = vci.VCIFile(filename_vci)

        if fasta_file.is_diploid():
            raise exceptions.G2GFastaError("Diploid Fasta files are not currently supported for patch")

        full_file = True

        if regions:
            full_file = False
            if len(regions) > 5:
                LOG.info("Regions: {} (showing 1st five)".format(", ".join(l for l in map(str, regions[:5]))))
            else:
                LOG.info("Regions: {}".format(", ".join(l for l in map(str, regions))))

        else:
            regions = []
            for chrom in fasta_file.references:
                regions.append(g2g.Region(chrom, 1, fasta_file.get_reference_length(chrom)))

        all_params = []
        for region in regions:
            LOG.debug('region={}'.format(region))
            LOG.debug('region.original_base={}'.format(region.original_base))
            params = FastaPatchParams()
            params.input_region = region
            params.input_file = filename_fasta
            params.temp_dir = temp_directory
            params.vci_file = filename_vci
            params.reverse = reverse
            params.offset = 0 if region.start <= 1 else region.start - 1

            if vci_file.is_haploid():
                LOG.debug("VCI File is haploid")
                params.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region), extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params.output_region = g2g.Region(region.seq_id, region.start, region.end)
                    params.output_header = fasta.FastaHeader(region.seq_id, "{}:{}-{}".format(region.seq_id, region.start, region.end))
                else:
                    params.output_region = g2g.Region("{}:{}-{}".format(region.seq_id, region.start, region.end), region.start, region.end)
                    params.output_header = fasta.FastaHeader("{}:{}-{}".format(region.seq_id, region.start, region.end), None)

                params.vci_query = "{}:{}-{}".format(region.seq_id, region.start, region.end)
                all_params.append(params)
            else:
                LOG.debug("VCI File is diploid")
                params.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region)+"_L", extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params.output_region = g2g.Region(region.seq_id+"_L", region.start, region.end)
                    params.output_header = fasta.FastaHeader(region.seq_id+"_L", "{}:{}-{}".format(region.seq_id, region.start, region.end))
                else:
                    params.output_region = g2g.Region("{}_L:{}-{}".format(region.seq_id, region.start, region.end), region.start, region.end)
                    params.output_header = fasta.FastaHeader("{}_L:{}-{}".format(region.seq_id, region.start, region.end), None)

                params.vci_query = "{}_L:{}-{}".format(region.seq_id, region.start, region.end)
                all_params.append(params)

                params_r = copy.deepcopy(params)
                params_r.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region)+"_R", extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params_r.output_region = g2g.Region(region.seq_id+"_R", region.start, region.end)
                    params_r.output_header = fasta.FastaHeader(region.seq_id+"_R", "{}:{}-{}".format(region.seq_id, region.start, region.end))
                else:
                    params_r.output_region = g2g.Region("{}_R:{}-{}".format(region.seq_id, region.start, region.end), region.start, region.end)
                    params_r.output_header = fasta.FastaHeader("{}_R:{}-{}".format(region.seq_id, region.start, region.end), None)

                params_r.vci_query = "{}_R:{}-{}".format(region.seq_id, region.start, region.end)
                all_params.append(params_r)

        args = zip(all_params)
        pool = multiprocessing.Pool(num_processes)
        results = pool.imap_unordered(wrapper, args)

        # parse results
        total = 0
        mode = 'wb'

        for c in results:
            if c is not None:
                total += c.count
                g2g_utils.concatenate_files([c.output_file], filename_output, False, mode)
                g2g_utils.delete_file(c.output_file)
                g2g_utils.delete_index_files(c.output_file)
                mode = 'ab'

        LOG.info("Patched {0:,} SNPs total".format(total))

        if dump_fasta:
            g2g_utils.dump_file_contents(filename_output)
        else:
            # move temp to final destination
            if bgzip:
                # remove the fai
                LOG.debug("removing the FAI index for {0}".format(filename_output))
                g2g_utils.delete_index_files(filename_output)

                LOG.info("Compressing and indexing...")

                if filename_output.lower().endswith((".fa", ".fasta")):
                    g2g_utils.bgzip_and_index_file(filename_output, "{0}.gz".format(filename_output), delete_original=True, file_format='fa')
                elif filename_output.lower().endswith(".gz"):
                    g2g_utils.bgzip_and_index_file(filename_output, filename_output, delete_original=True, file_format='fa')

    except KeyboardInterrupt:
        raise exceptions.KeyboardInterruptError()
    except Exception as e:
        raise exceptions.G2GError(str(e))
    finally:
        # clean up the temporary files
        g2g_utils.delete_dir(temp_directory)
        LOG.info("Patch complete: {0}".format(g2g_utils.format_time(start, time.time())))



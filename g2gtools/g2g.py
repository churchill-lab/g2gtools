# -*- coding: utf-8 -*-

from collections import OrderedDict
import os
import sys
import time

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from .bamsam import convert_bam_file
from .bed import convert_bed_file
from .chain import CHAIN_STRING, ChainFile, ChainIter, ChainEntry, collapse_entries
from .exceptions import G2GChainFileError, G2GVCFError, G2GValueError, G2GLocationError, G2GFastaError, G2GError
from .gtf import convert_gtf_file
from .gtf_db import get_genes_simple, get_transcripts_simple, gtf2db, get_genes, get_gene
from .g2g_utils import format_time, get_logger, parse_location, Location, merge_dicts, reverse_complement_sequence, wrap_sequence
import g2g_fileutils as g2g_fu
from .seq_offsets import offset_parse_chromosomes
from .vcf import VCF, parse_gt_new

from pysam import FastaFile
from pysam import TabixFile
import pysam

LOG = get_logger()


def file_convert(chain_file, input_file, output_file=None, file_format=None, reverse=False):
    """
    Convert a file

    :param chain_file:
    :param input_file:
    :param output_file:
    :param file_format:
    :return:
    """
    start = time.time()

    if file_format:
        file_format = file_format.upper()
        if file_format not in ['BED', 'BAM', 'SAM', 'GTF']:
            raise G2GValueError("Only BAM/SAM to BAM/SAM, GTF to GTF, or BED to BED are supported")
    else:
        # try to determine the type from the input
        file_all_caps = input_file.upper()
        if file_all_caps.endswith(('BAM', 'SAM')):
            file_format = 'BAM'
        elif file_all_caps.endswith('BED'):
            file_format = 'BED'
        elif file_all_caps.endswith('GTF'):
            file_format = 'GTF'
        else:
            raise G2GValueError("File format cannot be determined, please specify.")

    if file_format in ['BAM', 'SAM']:
        convert_bam_file(chain_file=chain_file, file_in=input_file, file_out=output_file, reverse=reverse)
    elif file_format in ['GTF']:
        convert_gtf_file(chain_file=chain_file, input_file=input_file, output_file=output_file, reverse=reverse)
    elif file_format in ['BED']:
        convert_bed_file(chain_file=chain_file, input_file=input_file, output_file=output_file, reverse=reverse)
    else:
        raise G2GValueError("Only BAM/SAM to BAM/SAM, GTF to GTF, or BED to BED are supported")

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def gtf2chain(chain_file, input_file, output_file, chain_genes=False):
    """

    :param chain_file:
    :param input_file:
    :param output_file:
    :param chain_genes:
    :return:
    """
    start = time.time()
    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))

    chain_file = g2g_fu.check_file(chain_file)
    input_file = g2g_fu.check_file(input_file)
    output_file = g2g_fu.check_file(output_file, 'w')
    output_file_dir = os.path.dirname(output_file)

    LOG.info("GTF FILE: {0}".format(input_file))
    LOG.info("FROM CHAIN FILE: {0}".format(chain_file))
    LOG.info("TO CHAIN FILE: {0}".format(output_file))

    temp_db = g2g_fu.gen_file_name("_g2gtempfile", output_file_dir, ".db3")

    gtf2db(input_file, temp_db)

    db2chain(chain_file, temp_db, output_file, chain_genes)

    g2g_fu.delete_file(temp_db)

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def convert_gene_transcripts(gene, chain_file, output_file):
    """
    Convert a gene's transcripts from the GTG DB file and chain file.

    :param gene: the gene information
    :param chain_file: the chain file
    :param output_file: the new chain file
    :return:
    """
    LOG.debug("Gene = {0}".format(str(gene)))

    if not isinstance(chain_file, ChainFile):
        chain_file = ChainFile(chain_file)

    transcripts = gene['transcripts']
    for t in transcripts:
        transcript = transcripts[t]
        LOG.debug("Transcript = {0}".format(str(transcript)))

        chain_entries = []
        from_start = None
        to_start = None

        for exon_id in transcript['exons']:
            exon = transcript['exons'][exon_id]
            LOG.debug("Exon = {0}".format(str(exon)))

            mappings = chain_file.find_mappings(gene['chrom'], exon['start'], exon['end'])

            if mappings and len(mappings) > 0:
                if not from_start:
                    from_start = mappings[0].from_start
                    to_start = mappings[0].to_start

                if len(mappings) == 1:
                    m = mappings[0]
                    c = ChainEntry()
                    #c.from_start = m.from_start
                    #c.from_end = m.from_end
                    #c.to_start = m.to_start
                    #c.to_end = m.to_end
                    c.lines.append([m.from_end - m.from_start])
                    chain_entries.append(c)
                else:
                    c = ChainEntry()
                    #c.from_start = mappings[0].from_start
                    #c.from_end = mappings[-1].from_end
                    #c.to_start = mappings[0].to_start
                    #c.to_end = mappings[-1].to_end

                    prev_mapping = None
                    sum_size = 0
                    sum_dq = 0
                    sum_dt = 0
                    dq = 0
                    prev_dq = 0
                    dt = 0
                    prev_dt = 0

                    for m in mappings:
                        if not prev_mapping:
                            prev_mapping = m
                        else:
                            prev_dt = dt
                            prev_dq = dq
                            chain_size = prev_mapping.from_end - prev_mapping.from_start
                            dt = m.from_start - prev_mapping.from_end
                            dq = m.to_start - prev_mapping.to_end

                            if dt > 0:
                                chain_size += prev_dq

                            sum_size += chain_size
                            sum_dq += dq
                            sum_dt += dt

                            c.lines.append([chain_size, dt, dq])
                            LOG.debug(c.lines[-1])
                            prev_mapping = m

                    chain_size = mappings[-1].from_end - mappings[-1].from_start
                    if dt > 0:
                        chain_size += dq
                        sum_size += dq

                    c.lines.append([chain_size])
                    chain_entries.append(c)

        # collapse exons
        if chain_entries and len(chain_entries) > 0:
            chain_entries = collapse_entries(chain_entries)
            sum_size = 0
            sum_dq = 0
            sum_dt = 0
            lines = []

            for line in chain_entries[0].lines:
                sum_size += line[0]
                if len(line) > 1:
                    sum_dq += line[1]
                    sum_dt += line[2]
                lines.append('\t'.join(map(str, line)))

            if output_file:
                outf = open(output_file, "a")
            else:
                outf = sys.stdout

            outf.write(CHAIN_STRING.format(CHAIN_STRING,
                        from_chr=t, from_length=sum_size + sum_dq,
                        from_start=0, from_end=sum_size + sum_dq,
                        to_chr=t, to_length=sum_size + sum_dt,
                        to_start=0, to_end=sum_size + sum_dt, id=t))
            outf.write("\n")
            outf.write("\n".join(lines))
            outf.write("\n")

            outf.close()


def db2chain(chain_file, input_file, output_file, chain_genes=False):
    """

    :param chain_file:
    :param input_file:
    :param output_file:
    :param chain_genes:
    :return:
    """
    start = time.time()

    if not isinstance(chain_file, ChainFile):
        chain_file = g2g_fu.check_file(chain_file)

    input_file = g2g_fu.check_file(input_file)
    output_file_name = g2g_fu.check_file(output_file, 'w')

    LOG.info("CHAIN FILE: {0}".format(chain_file))
    LOG.info("INPUT FILE: {0}".format(input_file))
    LOG.info("OUTPUT FILE: {0}".format(output_file_name))

    if chain_genes:
        LOG.info("CHAIN TYPE: GENES")
    else:
        LOG.info("CHAIN TYPE: TRANSCRIPTS")

    if not isinstance(chain_file, ChainFile):
        LOG.info("Parsing chain file...")
        chain_file = ChainFile(chain_file)
        LOG.info("Chain file parsed")

    LOG.info('Creating new chain file...')

    if chain_genes:
        LOG.debug("Generating chain for genes")

        for chromosome in chain_file.get_seqids():
            LOG.debug("Generating chain for genes in chromosome {0}".format(chromosome))

            for i, gene in enumerate(get_genes_simple(input_file, location=Location(chromosome))):
                LOG.debug("\n{0}".format(gene))
                chain_entries = []
                from_start = None
                to_start = None
                from_end = None
                to_end = None

                mappings = chain_file.find_mappings(gene.seqid, gene.start, gene.end)

                if gene.strand == 1:
                    if mappings and len(mappings) > 0:
                        if not from_start:
                            from_start = mappings[0].from_start
                            to_start = mappings[0].to_start

                        if len(mappings) == 1:
                            m = mappings[0]
                            c = ChainEntry()
                            c.lines.append([m.from_end - m.from_start])
                            chain_entries.append(c)
                        else:
                            c = ChainEntry()

                            prev_mapping = None
                            sum_size = 0
                            sum_dq = 0
                            sum_dt = 0
                            dq = 0
                            prev_dq = 0
                            dt = 0
                            prev_dt = 0

                            for m in mappings:
                                if not prev_mapping:
                                    prev_mapping = m
                                else:
                                    prev_dt = dt
                                    prev_dq = dq
                                    chain_size = prev_mapping.from_end - prev_mapping.from_start
                                    dt = m.from_start - prev_mapping.from_end
                                    dq = m.to_start - prev_mapping.to_end

                                    if dt > 0:
                                        chain_size += prev_dq

                                    sum_size += chain_size
                                    sum_dq += dq
                                    sum_dt += dt

                                    c.lines.append([chain_size, dt, dq])
                                    LOG.debug(c.lines[-1])
                                    prev_mapping = m

                            chain_size = mappings[-1].from_end - mappings[-1].from_start
                            if dt > 0:
                                chain_size += dq
                                sum_size += dq

                            c.lines.append([chain_size])
                            chain_entries.append(c)
                else:
                    if mappings and len(mappings) > 0:
                        if not from_end:
                            from_end = mappings[-1].from_end
                            to_end = mappings[-1].to_end

                        if len(mappings) == 1:
                            m = mappings[0]
                            c = ChainEntry()
                            c.lines.append([m.from_end - m.from_start])
                            chain_entries.append(c)
                        else:
                            c = ChainEntry()

                            prev_mapping = None
                            sum_size = 0
                            sum_dq = 0
                            sum_dt = 0
                            dq = 0
                            prev_dq = 0
                            dt = 0
                            prev_dt = 0

                            # reverse
                            mappings = mappings[::-1]

                            for m in mappings:
                                LOG.debug("CURRENT MAPPING: {0}".format(m))
                                if not prev_mapping:
                                    prev_mapping = m
                                else:
                                    LOG.debug("PREV MAPPING: {0}".format(prev_mapping))
                                    prev_dt = dt
                                    prev_dq = dq
                                    chain_size = prev_mapping.from_end - prev_mapping.from_start
                                    #dt = m.from_start - prev_mapping.from_end
                                    #dq = m.to_start - prev_mapping.to_end
                                    dt = prev_mapping.from_start - m.from_end
                                    dq = prev_mapping.to_start - m.to_end
                                    LOG.debug("dt={0}, dq={1}".format(dt, dq))

                               #if dt > 0:
                                    #    LOG.debug("DT > 0, ADDING to current chain_size {0}".format(chain_size))
                                    #    chain_size += prev_dq

                                    sum_size += chain_size
                                    sum_dq += dq
                                    sum_dt += dt

                                    c.lines.append([chain_size, dt, dq])
                                    LOG.debug(c.lines[-1])
                                    prev_mapping = m

                            LOG.debug("finding last...{0}".format(mappings[-1]))
                            chain_size = mappings[-1].from_end - mappings[-1].from_start
                            #if dt > 0:
                            #    LOG.debug("WHOA {0}".format(dt))
                            #    LOG.debug("DT > 0, ADDING to current chain_size {0}".format(chain_size))
                            #    chain_size += dq
                            #    sum_size += dq

                            c.lines.append([chain_size])
                            LOG.debug(c.lines[-1])
                            chain_entries.append(c)

                if chain_entries and len(chain_entries) > 0:
                    sum_size = 0
                    sum_dq = 0
                    sum_dt = 0
                    lines = []

                    for line in chain_entries[0].lines:
                        sum_size += line[0]
                        if len(line) > 1:
                            sum_dq += line[1]
                            sum_dt += line[2]
                        lines.append('\t'.join(map(str, line)))

                    if output_file:
                        outf = open(output_file, "a")
                    else:
                        outf = sys.stdout

                outf.write(CHAIN_STRING.format(CHAIN_STRING,
                            from_chr=gene.seqid, from_length=sum_size + sum_dq,
                            from_start=0, from_end=sum_size + sum_dq,
                            to_chr=gene.seqid, to_length=sum_size + sum_dt,
                            to_start=0, to_end=sum_size + sum_dt, id=gene.ensembl_id))
                outf.write("\n")
                outf.write("\n".join(lines))
                outf.write("\n")

                outf.close()
    else:
        for chromosome in chain_file.get_seqids():
            LOG.debug("Generating chain for transcripts in chromosome {0}".format(chromosome))

            for i, transcript in enumerate(get_transcripts_simple(input_file, location=Location(chromosome))):
                LOG.debug("Transcript = {0}".format(transcript))
                chain_entries = []
                from_start = None
                to_start = None
                from_end = None
                to_end = None
                transcript.exons = OrderedDict(sorted(transcript.exons.items(), key=lambda x: x[1].exon_number))

                for ensembl_id, exon in transcript.exons.iteritems():
                    LOG.debug("Exon = {0}".format(exon))

                    mappings = chain_file.find_mappings(exon.seqid, exon.start, exon.end)

                    if exon.strand == 1:
                        if mappings and len(mappings) > 0:
                            if not from_start:
                                from_start = mappings[0].from_start
                                to_start = mappings[0].to_start

                            if len(mappings) == 1:
                                m = mappings[0]
                                c = ChainEntry()
                                c.lines.append([m.from_end - m.from_start])
                                chain_entries.append(c)
                            else:
                                c = ChainEntry()

                                prev_mapping = None
                                sum_size = 0
                                sum_dq = 0
                                sum_dt = 0
                                dq = 0
                                prev_dq = 0
                                dt = 0
                                prev_dt = 0

                                for m in mappings:
                                    if not prev_mapping:
                                        prev_mapping = m
                                    else:
                                        prev_dt = dt
                                        prev_dq = dq
                                        chain_size = prev_mapping.from_end - prev_mapping.from_start
                                        dt = m.from_start - prev_mapping.from_end
                                        dq = m.to_start - prev_mapping.to_end

                                        if dt > 0:
                                            chain_size += prev_dq

                                        sum_size += chain_size
                                        sum_dq += dq
                                        sum_dt += dt

                                        c.lines.append([chain_size, dt, dq])
                                        LOG.debug(c.lines[-1])
                                        prev_mapping = m

                                chain_size = mappings[-1].from_end - mappings[-1].from_start
                                if dt > 0:
                                    chain_size += dq
                                    sum_size += dq

                                c.lines.append([chain_size])
                                chain_entries.append(c)
                    else:
                        if mappings and len(mappings) > 0:
                            if not from_end:
                                from_end = mappings[-1].from_end
                                to_end = mappings[-1].to_end

                            if len(mappings) == 1:
                                m = mappings[0]
                                c = ChainEntry()
                                c.lines.append([m.from_end - m.from_start])
                                chain_entries.append(c)
                            else:
                                c = ChainEntry()

                                prev_mapping = None
                                sum_size = 0
                                sum_dq = 0
                                sum_dt = 0
                                dq = 0
                                prev_dq = 0
                                dt = 0
                                prev_dt = 0

                                # reverse
                                mappings = mappings[::-1]

                                for m in mappings:
                                    LOG.debug("CURRENT MAPPING: {0}".format(m))
                                    if not prev_mapping:
                                        prev_mapping = m
                                    else:
                                        LOG.debug("PREV MAPPING: {0}".format(prev_mapping))
                                        prev_dt = dt
                                        prev_dq = dq
                                        chain_size = prev_mapping.from_end - prev_mapping.from_start
                                        #dt = m.from_start - prev_mapping.from_end
                                        #dq = m.to_start - prev_mapping.to_end
                                        dt = prev_mapping.from_start - m.from_end
                                        dq = prev_mapping.to_start - m.to_end
                                        LOG.debug("dt={0}, dq={1}".format(dt, dq))

                                   #if dt > 0:
                                        #    LOG.debug("DT > 0, ADDING to current chain_size {0}".format(chain_size))
                                        #    chain_size += prev_dq

                                        sum_size += chain_size
                                        sum_dq += dq
                                        sum_dt += dt

                                        c.lines.append([chain_size, dt, dq])
                                        LOG.debug(c.lines[-1])
                                        prev_mapping = m

                                LOG.debug("finding last...{0}".format(mappings[-1]))
                                chain_size = mappings[-1].from_end - mappings[-1].from_start
                                #if dt > 0:
                                #    LOG.debug("WHOA {0}".format(dt))
                                #    chain_size += dq
                                #    sum_size += dq

                                c.lines.append([chain_size])
                                LOG.debug(c.lines[-1])
                                chain_entries.append(c)

                # collapse exons
                if chain_entries and len(chain_entries) > 0:
                    LOG.debug('>>>>>>>')
                    for c in chain_entries:
                        LOG.debug(str(c))
                    LOG.debug('>>>>>>>')
                    chain_entries = collapse_entries(chain_entries)
                    sum_size = 0
                    sum_dq = 0
                    sum_dt = 0
                    lines = []

                    for line in chain_entries[0].lines:
                        sum_size += line[0]
                        if len(line) > 1:
                            sum_dq += line[1]
                            sum_dt += line[2]
                        lines.append('\t'.join(map(str, line)))

                    if output_file:
                        outf = open(output_file, "a")
                    else:
                        outf = sys.stdout

                    outf.write(CHAIN_STRING.format(CHAIN_STRING,
                                from_chr=transcript.seqid, from_length=sum_size + sum_dq,
                                from_start=0, from_end=sum_size + sum_dq,
                                to_chr=transcript.seqid, to_length=sum_size + sum_dt,
                                to_start=0, to_end=sum_size + sum_dt, id=transcript.ensembl_id))
                    outf.write("\n")
                    outf.write("\n".join(lines))
                    outf.write("\n")

                    outf.close()

    LOG.info('New chain file created')

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def offset_chromosome_to_chain(chromosome, chain_file):
    """
    Convert Seqnature offset file to chain file

    :param chromosome: the chromosome data
    :param chain_file: the chain file to write to
    """
    if not os.path.exists(chromosome['file_path']):
        raise IOError('Cannot find {0}'.format(chromosome['file_path']))

    prev_offset_pos = 0
    prev_offset_amt = 0
    first_offset = 0
    ref_after_pos = 0

    sum_size = 0
    sum_dt = 0
    sum_dq = 0
    chain_data = []

    fd = open(chromosome['file_path'], 'r')

    for line in fd:
        elem = line.strip().split()
        offset_pos = int(elem[0])
        offset_amt = int(elem[1])

        LOG.debug('LINE: {0}'.format(elem))

        chain_size = 0
        dt = 0
        dq = 0

        if prev_offset_pos == 0 and prev_offset_amt == 0:
            # first time through
            LOG.debug("ref_after_pos=" + str(ref_after_pos))
            LOG.debug('offset_pos=' + str(offset_pos) + '\toffset_amt=' + str(offset_amt))
            LOG.debug('prev_pos=' + str(prev_offset_pos) + '\tprev_amt=' + str(prev_offset_amt))
            offset_diff = offset_amt
            LOG.debug('offset_diff= ' + str(offset_diff))

            if offset_amt >= 0:
                dt = offset_amt
                ref_after_pos = offset_pos
            else:
                dq = offset_amt * -1
                ref_after_pos = offset_pos + dq

            chain_size = offset_pos - 1
            first_offset = chain_size
            sum_size = chain_size
            sum_dt = dt
            sum_dq = dq
            chain_data.append([chain_size, dt, dq])

            prev_offset_pos = offset_pos
            prev_offset_amt = offset_amt
            LOG.debug('>> {0}\t{1}\t{2}'.format(sum_size, sum_dt, sum_dq))
        else:
            offset_diff = offset_amt - prev_offset_amt

            LOG.debug("PREV ref_after_pos=" + str(ref_after_pos))
            LOG.debug("offset_pos=" + str(offset_pos))

            if offset_diff >= 0:
                dt = offset_diff
                chain_size = offset_pos - ref_after_pos
                ref_after_pos = offset_pos
            else:
                dq = offset_diff * -1
                chain_size = offset_pos - ref_after_pos
                ref_after_pos = offset_pos + dq

            LOG.debug("ref_after_pos=" + str(ref_after_pos))


            sum_size += chain_size
            sum_dt += dt
            sum_dq += dq

            # fix any 0 length
            if chain_size != 0:

                chain_data.append([chain_size, dt, dq])
                LOG.debug(chain_data[-1])
            else:
                LOG.debug("offset_diff=" + str(offset_diff))
                LOG.debug("offset_amt=" + str(offset_amt))
                LOG.debug("prev_offset_amt=" + str(prev_offset_amt))
                LOG.debug("COMBINING!!!!!!!!!!")
                LOG.debug(chain_data[-1])
                LOG.debug("dt={0}".format(dt))
                LOG.debug("dq={0}".format(dq))
                chain_data[-1][1] += dt
                chain_data[-1][2] += dq

            prev_offset_pos = offset_pos
            prev_offset_amt = offset_amt
            LOG.debug('>>SUM {0}\t{1}\t{2}'.format(sum_size, sum_dt, sum_dq))
            LOG.debug('>>    {0}\t{1}\t{2}'.format(chain_size, dt, dq))

    # not sure if the following assumption will work in all cases
    # the last line of the chain segment is the last fragment length
    if chromosome['from_length'] > chromosome['to_length']:
        last_chain = chromosome['from_length'] - (sum_size + sum_dt)
        sum_size += last_chain
    else:
        last_chain = chromosome['to_length'] - (sum_size + sum_dq)
        sum_size += last_chain
    chain_fd = sys.stdout
    if chain_file:
        chain_fd = open(chain_file, "a")

    chain_fd.write(CHAIN_STRING.format(CHAIN_STRING,
            from_chr=chromosome['from_chr'], from_length=chromosome['from_length'],
            from_start=0, from_end=sum_size + sum_dt,
            to_chr=chromosome['to_chr'], to_length=chromosome['to_length'],
            to_start=0, to_end=sum_size + sum_dq, id=chromosome['to_chr']))
    chain_fd.write("\n")

    for c in chain_data:
        chain_fd.write("{0}\t{1}\t{2}\n".format(c[0], c[1], c[2]))

    chain_fd.write("{0}\n\n".format(last_chain))

    chain_fd.close()


def offset2chain(from_file, to_file, output_file):
    """
    Convert Seqnature offset files to new chain file.

    :param from_file: from Chromosome File (see docs)
    :param to_file: to Chromosome File (see docs)
    :param output_file: the output chain file
    """
    start = time.time()

    from_file = g2g_fu.check_file(from_file)
    to_file = g2g_fu.check_file(to_file)

    output_file_name = g2g_fu.check_file(output_file, 'w')
    g2g_fu.delete_file(output_file_name)

    LOG.info("FROM FILE: {0}".format(from_file))
    LOG.info("TO FILE: {0}".format(to_file))
    LOG.info("CHAIN FILE: {0}".format(output_file_name))

    LOG.info("Generating chain file...")

    try:
        chromosomes = offset_parse_chromosomes(from_file, to_file)

        for c, chromosome in chromosomes.iteritems():
            LOG.debug('Examining chromosome: {0}'.format(chromosome))
            if chromosome['file_path']:
                offset_chromosome_to_chain(chromosome, output_file)
            else:
                LOG.debug("No file for {0}, so skipping".format(chromosome))

        LOG.info("Chain file created")

    except Exception, e:
        raise G2GChainFileError("Unable to generate chain file")

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def fasta_extract_transcripts(fasta_file, database_file, output, raw=False):
    start = time.time()

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_fu.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file = g2g_fu.check_file(database_file)

    fasta_out = sys.stdout

    if output:
        output = g2g_fu.check_file(output, 'w')
        fasta_out = open(output, "w")

    LOG.info("FASTA FILE: {0}".format(fasta.filename))
    LOG.info("DATABASE FILE: {0}".format(database_file))
    LOG.info("OUTPUT FILE: {0}".format(fasta_out.name))

    try:
        transcripts = get_transcripts_simple(database_file)

        for i, transcript in enumerate(transcripts):
            LOG.debug("Transcript={0}".format(transcript))

            if transcript.seqid not in fasta.references:
                continue

            new_sequence = StringIO()
            locations = []
            for ensembl_id, exon in transcript.exons.iteritems():
                LOG.debug("Exon ID={0};{1}".format(ensembl_id, exon))

                partial_seq = fasta.fetch(exon.seqid, exon.start-1, exon.end)
                partial_seq_str = str(partial_seq)

                if transcript.strand == 1:
                    new_sequence.write(partial_seq_str)
                else:
                    partial_seq_str = str(reverse_complement_sequence(partial_seq))
                    new_sequence.write(partial_seq_str)

                LOG.debug("{0}:{1}-{2} (Length: {3})\n{4}".format(exon.seqid, exon.start, exon.end, len(partial_seq), partial_seq_str))
                locations.append("{0}:{1}-{2}".format(exon.seqid, exon.start, exon.end))

            if raw:
                fasta_out.write(new_sequence.getvalue())
            else:
                fasta_id = ">{0} {1}|{2}\n".format(transcript.ensembl_id, '-' if transcript.strand == -1 else '+', "|".join(locations))
                fasta_out.write(fasta_id)

                for line in wrap_sequence(new_sequence.getvalue()):
                    fasta_out.write(line.strip())
                    fasta_out.write('\n')

    except G2GValueError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def fasta_extract_exons(fasta_file, database_file, output, raw=False):
    start = time.time()

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_fu.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file = g2g_fu.check_file(database_file)

    fasta_out = sys.stdout

    if output:
        output = g2g_fu.check_file(output, 'w')
        fasta_out = open(output, "w")

    LOG.info("FASTA FILE: {0}".format(fasta.filename))
    LOG.info("DATABASE FILE: {0}".format(database_file))
    LOG.info("OUTPUT FILE: {0}".format(fasta_out.name))

    try:
        transcripts = get_transcripts_simple(database_file)
        for i, transcript in enumerate(transcripts):

            if transcript.seqid not in fasta.references:
                continue

            for ensembl_id, exon in transcript.exons.iteritems():
                LOG.debug("Exon={0}".format(exon))

                partial_seq = fasta.fetch(exon.seqid, exon.start-1, exon.end)
                partial_seq_str = partial_seq

                if transcript.strand == -1:
                    partial_seq_str = str(reverse_complement_sequence(partial_seq))

                LOG.debug("{0}:{1}-{2} (Length: {3})\n{4}".format(exon.seqid, exon.start, exon.end, len(partial_seq), partial_seq_str))

                if raw:
                    fasta_out.write(partial_seq_str)
                else:
                    fasta_id = ">{0} {1}:{2}-{3}\n".format(exon.ensembl_id, exon.seqid, exon.start, exon.end)
                    fasta_out.write(fasta_id)

                    for line in wrap_sequence(partial_seq_str):
                        fasta_out.write(line.strip())
                        fasta_out.write('\n')

    except G2GValueError as e:
        LOG.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        LOG.info(e.msg.rstrip())
        raise e

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def dbfetch(input_file, output_file=None, locations=None, ids=None, display_genes=True, display_transcripts=True, display_exons=True, overlap=True):
    """
    """
    start = time.time()

    if input_file and os.path.exists(input_file):
        input_file = os.path.abspath(input_file)
    else:
        raise G2GValueError("The following G2G DB file does not exist: {0}".format(input_file))

    output_file_name = None

    if output_file:
        output_file = os.path.abspath(output_file)
        output_file_dir = os.path.dirname(output_file)
        output_file_name = output_file

        if not os.access(output_file_dir, os.W_OK | os.X_OK):
            raise G2GValueError("Cannot generate file in specified directory: {0}".format(output_file_dir))

    LOG.info("DB FILE: {0}".format(input_file))

    if output_file_name:
        LOG.info("OUTPUT FILE: {0}".format(output_file_name))

    if (locations is None and ids is None) or (locations and ids):
        raise G2GValueError("Can only search by location or ids separately")

    genes = {}
    if locations:
        try:
            for l in locations:
                if not isinstance(l, Location):
                    l = parse_location(l)

                genes = merge_dicts(genes, get_genes(input_file, l))
        except G2GLocationError, le:
            LOG.debug("Unable to parse location, {0}".format(le.message))
            raise le
    else:
        for i in ids:
            genes = merge_dicts(genes, get_gene(input_file, i))

    for k,v in genes.iteritems():
        indent = 0
        if display_genes:
            print v

        for k1, v1 in v.transcripts.iteritems():
            if display_transcripts:
                indent = 1 if display_genes else 0
                print '\t'*indent, str(v1)

            if display_exons:
                if display_genes and display_transcripts:
                    indent = 2
                elif display_genes or display_transcripts:
                    indent = 1
                else:
                    indent = 0

                for k2,v2 in v1.exons.iteritems():
                    print '\t'*indent, str(v2)

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


def fasta_transform(fasta_file, chain_file, locations, output_file, bgzip=False, reverse=False):
    """

    :param fasta_file:
    :param chain_file:
    :param locations:
    :param output_file:
    :param bgzip:
    :param reverse:
    :return:
    """
    start = time.time()

    if not isinstance(fasta_file, FastaFile):
        fasta_file = g2g_fu.check_file(fasta_file)

    if not isinstance(chain_file, ChainIter):
        chain_file = g2g_fu.check_file(chain_file)

    output_file = g2g_fu.check_file(output_file, 'w')
    g2g_fu.delete_file(output_file)
    g2g_fu.delete_index_files(output_file)

    LOG.info("FASTA FILE: {0}".format(fasta_file))
    LOG.info("CHAIN FILE: {0}".format(chain_file))
    LOG.info("OUTPUT FILE: {0}".format(output_file))
    LOG.info("BGZIP: {0}".format(bgzip))
    LOG.info("REVERSE: {0}".format(reverse))

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta = FastaFile(fasta_file)

    if not isinstance(chain_file, ChainIter):
        chain_file = ChainIter(chain_file, reverse=reverse)

    seq_ids = []

    if locations:
        LOG.debug("Have locations")
        new_locations = []
        for l in locations:
            if isinstance(l, Location):
                new_locations.append(l)
            else:
                new_locations.append(parse_location(l))
            seq_ids.append(new_locations[-1].seqid)
        locations = new_locations
    else:
        LOG.debug("Calculating locations")
        locations = [parse_location("{0}:1-{1}".format(a, fasta.get_reference_length(a)), 1) for a in fasta.references]
        seq_ids = [a for a in fasta.references]

    temp_output_file = output_file

    if bgzip:
        if g2g_fu.get_extension(output_file) != 'gz':
            output_file = "{0}.gz".format(output_file)
        else:
            temp_output_file = temp_output_file[:-3]

    fasta_out = open(temp_output_file, "w")

    LOG.info("Transforming...")

    chr_info = {}

    try:
        # will need a better way, but for now...
        LOG.info("Parsing chain file...")
        for line in chain_file:
            if len(line) > 7:
                LOG.debug("Adding chromosome {0}".format(chain_file.current_chain_header[1]))
                chr_info[chain_file.current_chain_header[1]] = {'from_size': line[2], 'from_start': line[4], 'from_end': line[5],
                                  'to_size': line[7], 'to_start': line[9], 'to_end': line[10],
                                  'header_chain':chain_file.current_chain_header, 'lines': []}
            else:
                chr_info[chain_file.current_chain_header[1]]['lines'].append(line)

        LOG.info("Chain file parsed")

        insertion_bases = 0
        deletion_bases = 0

        for location in locations:
            LOG.info("Processing chromosome={0}".format(location.seqid))
            LOG.debug("Location: {0}".format(location))

            chrom_size_from = chr_info[location.seqid]['from_size']
            chrom_size_to = chr_info[location.seqid]['to_size']

            last_pos = chr_info[location.seqid]['from_start']
            new_sequence = StringIO()
            chain_file.reset()

            for chain_line in chr_info[location.seqid]['lines']:
                LOG.debug("\nLINE: {0} : {1}".format(chain_file.line_no, chain_line))

                if len(chain_line) == 1:
                    # last line
                    fragment = chain_line[0]

                    partial_seq = fasta.fetch(location.seqid, last_pos, last_pos + fragment)
                    new_sequence.write(str(partial_seq))

                    if len(new_sequence.getvalue()) < chrom_size_to:
                        LOG.warn("Length's do not match, chromosome length in chain: {0}, sequence length: {1}".format(chrom_size_to, len(new_sequence.getvalue())))

                    fasta_out.write(">{0} {1}:{2}-{3}\n".format(location.seqid, location.seqid, chr_info[location.seqid]['from_start'] + 1, chrom_size_to))

                    for l in wrap_sequence(new_sequence.getvalue()):
                        fasta_out.write(l.strip())
                        fasta_out.write('\n')

                    break

                else:

                    # fragment_size dt_size dq_size same_bases dt_bases dq_bases

                    fragment = chain_line[0]
                    dt = chain_line[1 if not reverse else 2]
                    dq = chain_line[2 if not reverse else 1]
                    same = chain_line[3]
                    dt_bases = chain_line[4 if not reverse else 5]
                    dq_bases = chain_line[5 if not reverse else 4]

                    partial_seq = fasta.fetch(location.seqid, last_pos, last_pos + fragment)
                    new_sequence.write(partial_seq)

                    if dq > 0:
                        # insertion
                        LOG.debug("INSERTION")
                        new_sequence.write(dq_bases)
                        LOG.debug("{0}:{1}-{2} (Length: {3})".format(location.seqid, last_pos, last_pos + fragment, len(partial_seq)))
                        if len(partial_seq) > 100:
                            LOG.debug("{0}...{1}".format(partial_seq[:10], partial_seq[-10:]))
                        else:
                            LOG.debug(partial_seq)
                        LOG.debug("Adding {0}".format(dq_bases))
                        LOG.debug("SAME={0}, {1}".format(same, partial_seq[-(len(same)):]))

                        insertion_bases += dq

                    if dt > 0:
                        # deletion
                        LOG.debug("DELETION")
                        last_pos += dt
                        LOG.debug("skipping ahead {0} bases".format(dt))

                        deletion_bases += dt

                    last_pos += fragment

                    LOG.debug("LAST_POS={0}, INSERTIONS={1}, DELETIONS={2}, DIFF={3}".format(last_pos, insertion_bases, deletion_bases, (insertion_bases - deletion_bases)))

        # bgzip and index
        if bgzip:
            LOG.info("Compressing and indexing...")
            g2g_fu.bgzip_index(temp_output_file, output_file, 'fa')

    except G2GLocationError, le:
        LOG.debug("Unable to parse location, {0}".format(le.message))
        raise le
    except G2GValueError as e:
        LOG.debug("Unable to parse alocation, {0}".format(e.message))
        raise e
    except G2GFastaError as e:
        LOG.debug("Unable to parse blocation, {0}".format(e.message))
        raise e
    except TypeError as e:
        LOG.debug("Unable to parse clocation, {0}".format(e.message))
        raise e
        #print str(e)
        #raise G2GError("Improper chain file")

    LOG.info("Transforming complete")

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))

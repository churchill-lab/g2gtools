# -*- coding: utf-8 -*-

import argparse
import sys
import time

from . import bsam
from . import bed
from . import exceptions
from . import fasta
from . import fasta_patch
from . import fasta_transform
from . import g2g
from . import g2g_utils
from . import gtf
from . import gtf_db
from . import vci
from . import vcf2vci

LOG = g2g.get_logger()

#TODO: Ask KB about --pass, --quality  do we want here and also at the transform and patch level?

#TODO: Assumption.., extract and location work on INPUT genome, but shouldn't we be able to say, I know what the position is in CAST of where I want to look... show me this sequence in the ref


def command_convert(raw_args, prog=None):
    """
    Convert a file

    Usage: convert [-options] -v <VCI file> -i <input file>

    Required Parameters:
    -i, --input <input file>         input file to convert
    -v, --vci <VCI file>             input VCI file used in conversion

    Optional Parameters:
    -f, --format <bam|sam|gtf|bed>   file format
    -r, --reverse                    reverse the Chain file
    -o, --output <output file>       file to create

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    If no output file is specified [-o], the converted information will be redirected
    to standard out.

    To specify the type of file, use the [-t] option.

    Only the conversion from the same type of file is supported.

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_convert.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-v", "--vci", dest="vci", metavar="VCI_File")

    # optional
    parser.add_argument("-f", "--format", dest="format", metavar='bam|sam|gtf|bed', choices=['bam', 'sam', 'gtf', 'bed'])
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")
    parser.add_argument("-r", "--reverse", dest="reverse", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if args.help:
        g2g.exit("", parser)

    if not args.input:
        g2g.exit("No input file was specified.", parser)

    if not args.vci:
        g2g.exit("No VCI file was specified.", parser)

    try:
        file_format = None

        if args.format:
            file_format = file_format.upper()
            if file_format not in ['BED', 'BAM', 'SAM', 'GTF']:
                raise exceptions.G2GValueError("Only BAM/SAM to BAM/SAM, GTF to GTF, or BED to BED are supported")
        else:
            # try to determine the type from the input
            file_all_caps = args.input.upper()
            if file_all_caps.endswith(('BAM', 'SAM')):
                file_format = 'BAM'
            elif file_all_caps.endswith('BED'):
                file_format = 'BED'
            elif file_all_caps.endswith('GTF'):
                file_format = 'GTF'
            else:
                raise exceptions.G2GValueError("File format cannot be determined, please specify.")

        if file_format in ['BAM', 'SAM']:
            bsam.convert_bam_file(vci_file=args.vci, file_in=args.input, file_out=args.output, reverse=args.reverse)
        elif file_format in ['GTF']:
            gtf.convert_gtf_file(vci_file=args.vci, input_file=args.input, output_file=args.output, reverse=args.reverse)
        elif file_format in ['BED']:
            bed.convert_bed_file(vci_file=args.vci, input_file=args.input, output_file=args.output, reverse=args.reverse)
        else:
            raise exceptions.G2GValueError("Only BAM/SAM to BAM/SAM, GTF to GTF, or BED to BED are supported")
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as ve:
        g2g.exit(ve, parser)
    except exceptions.G2GChainFileError as cfe:
        g2g.exit(cfe, parser)
    except exceptions.G2GBAMError as e:
        g2g.exit(e, parser)
    except exceptions.G2GBedError as e:
        g2g.exit(e, parser)


def command_vcf2vci(raw_args, prog=None):
    """
    Convert snp and indel VCF file(s) to a VCI file

    Usage: vcf2vci [-options] [-v <indel VCF file> *] -s <strain> -o <output VCI file>

    Required Parameters:
        -o, --output <output file>       VCI file to create
        -s, --strain <strain>            name of strain column in VCF file
        -v, --vcf <vcf_file>             snp/indel VCF file

    Optional Parameters:
        --nobgzip                        DO NOT compress and index output
        --diploid                        ignore hetorzygotes
        --keep                           keep track of VCF lines that cannot be converted to Chain file
        --pass                           use only VCF lines that have a PASS for the filter value
        --quality                        filter on quality, FI=PASS
        -n, --num_processes <number>     the number of processes to use, defaults to the number of cores

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_vcf2vci.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-o", "--output", dest="output", metavar="VCI_File")
    parser.add_argument("-s", "--strain", dest="strain", metavar="strain")
    parser.add_argument("-v", "--vcf", dest="vcf_files", metavar="vcf_file", action='append')

    # optional
    parser.add_argument("--nobgzip", dest="nobgzip", action='store_false', default=True)
    parser.add_argument("--diploid", dest="diploid", action='store_true')
    parser.add_argument("--keep", dest="keep", action='store_true')
    parser.add_argument("--pass", dest="passed", action='store_true')
    parser.add_argument("--quality", dest="quality", action='store_true')
    parser.add_argument("-n", "--numprocesses", type=int, dest="numprocesses", metavar="number_of_processes")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    if not args.vcf_files:
        g2g.exit("No VCF file was specified.", parser)

    if not args.output:
        g2g.exit("No output VCF file was specified.", parser)

    if not args.strain:
        g2g.exit("No strain was specified.", parser)

    try:
        vcf2vci.process(args.vcf_files, args.output, args.strain, args.keep, args.passed, args.quality, args.diploid, args.numprocesses, args.nobgzip)
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as e:
        g2g.exit(e, parser)
    except exceptions.G2GVCFError as e:
        g2g.exit(e, parser)


def command_fasta_extract(raw_args, prog=None):
    """
    Extract a sequence from a Fasta file give a location.

    Usage: extract [-options] -i <Fasta file> [-l <location> | -b <BED file> | -db <Database> | -id <id>]

    Required Parameters:
        -i, --input <Fasta file>         Fasta file to extract sequence from, BGZIP files supported

    Optional Parameters:
        -b, --bed <BED file>             BED file
            - OR -
        -id, --identifier <identifier>   Fasta identifier
            - OR -
        -r, --region <region>            seqid:start:end

        -db, --database <DB file>        Database file
            --exons                      For use with -db, defaults to transcript extraction
            --genes                      For use with -db, defaults to transcript extraction
            --transcripts                For use with -db, defaults to transcript extraction

        -v, --vci <VCI file>             input VCI file, matching input Fasta file

        --complement                     return complement of sequence
        --reverse                        return reverse of sequence
        --reversecomplement              return reverse complement of sequence

        --raw                            show just sequence

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    Note: Locations specified on command line are 1 based coordinates.
          Locations specified via BED file are 0 based coordinates.

    """
    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_fasta_extract.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")

    # optional
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-id", "--identifier", dest="id", metavar="identifier")
    parser.add_argument("-r", "--region", dest="region", metavar="chr1:100000-200000")

    parser.add_argument("-db", "--database", dest="database", metavar="Database")
    parser.add_argument("--exons", dest="exons", action="store_true", default=False)
    parser.add_argument("--genes", dest="genes", action="store_true", default=False)
    parser.add_argument("--transcripts", dest="transcripts", action="store_true", default=False)

    parser.add_argument("--complement", dest="complement", action="store_true", default=False)
    parser.add_argument("--reverse", dest="reverse", action="store_true", default=False)
    parser.add_argument("--reversecomplement", dest="reversecomplement", action="store_true", default=False)

    parser.add_argument("-v", "--vci", dest="vci", metavar="VCI_File")
    parser.add_argument("-vr", "--vci-reverse", dest="vci_reverse", action="store_true", default=False)

    parser.add_argument("--raw", dest="raw", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    if not args.fasta:
        g2g.exit("No Fasta file was specified.", parser)

    flag_count = 0
    flag_count += (1 if args.complement else 0)
    flag_count += (1 if args.reverse else 0)
    flag_count += (1 if args.reversecomplement else 0)

    if flag_count > 1:
        g2g.exit("Please specify only one of: --complement, --reverse, --reversecomplement", parser)

    reverse = (args.reverse or args.reversecomplement)
    complement = (args.complement or args.reversecomplement)

    if args.database:
        # allow only either exon, gene, or transcript extraction, bot all
        output_count = 0
        output_count += (1 if args.genes else 0)
        output_count += (1 if args.transcripts else 0)
        output_count += (1 if args.exons else 0)

        if output_count != 1:
            g2g.exit("Please specify one of: --exons, --genses, or --transcripts", parser)

    location_count = 0
    location_count += (1 if args.region else 0)
    location_count += (1 if args.bed_file else 0)
    location_count += (1 if args.database else 0)
    location_count += (1 if args.id else 0)

    if location_count > 1:
        g2g.exit("Please specify only one of: -r, -b, -id, or -db", parser)



    all_regions = []
    transcript_info = None

    # parse the regions we are extracting

    try:
        fasta_file = g2g_utils.check_file(args.fasta)
        fasta_file = fasta.FastaFile(fasta_file)

        # get the regions if there are some

        if args.region:
            # simple region
            LOG.debug('ARGS.REGION={}'.format(args.region))
            region = g2g.parse_region(args.region, base=1)


            #region.name = "{}:{}-{}".format(region.seq_id, region.start, region.end)

            LOG.debug("--> start = {}".format(region.start))
            LOG.debug("--> _start = {}".format(region._start))
            LOG.debug("--> get_start() = {}".format(region.get_start()))

            LOG.debug('REGION PARSED={}'.format(region))
            all_regions.append(region)

            # TEMPORARY OVERRIDE
            #all_regions = [g2g.parse_region(args.region)]

        elif args.bed_file:
            # bed file
            bed_file = bed.BED(args.bed_file)
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else '+'
                    LOG.debug(bed_rec)
                    all_regions.append(g2g.Region(bed_rec.chrom, bed_rec.start, bed_rec.end, strand, bed_rec.name, 0))

        # on the fly
        if args.vci:

            # let's eliminate what isn't supported

            if fasta_file.is_diploid():
                g2g.exit("Diploid Fasta files are not supported with -v, --vci", parser)
            elif args.id:
                g2g.exit("Options -v, --vci cannot be used with --id, --identifier", parser)
            elif args.database:
                g2g.exit("Options -v, --vci cannot be used with --db, --database", parser)
            elif flag_count >= 1:
                g2g.exit("Options --complement, --reverse, --reversecomplement cannot be used with VCI.", parser)

            #fasta_onthefly.process(filename_fasta=args.fasta, filename_vci=args.vci, regions=all_regions,
            #                        filename_output=None, bgzip=False, reverse=args.vci_reverse, num_processes=None,
            #                        transcript_info=transcript_info)

            if args.region:
                #overridden here, fasta.extract expects 1 base

                all_regions = [g2g.parse_region(args.region)]

            fasta_transform.process(filename_fasta=args.fasta, filename_vci=args.vci, regions=all_regions,
                                    filename_output=None, bgzip=False, reverse=args.vci_reverse, num_processes=None, also_patch=True)

        # normal extraction
        else:
            if args.id:
                fasta.extract_id(args.fasta, args.id, output=None, reverse=reverse, complement=complement, raw=args.raw)
                return

            elif args.database:
                if flag_count >= 1:
                    g2g.exit("Options --complement, --reverse, --reversecomplement cannot be used with Database.", parser)

                if args.transcripts:
                    fasta.fasta_extract_transcripts(fasta_file, args.database, None, raw=args.raw)
                    return

                elif args.genes:
                    LOG.info("Extracting genes from database...")
                    genes = gtf_db.get_genes_simple(args.database)

                    LOG.debug("Genes extracted, transforming to locations")
                    for gene in genes:
                        r = g2g.Region(gene.seqid, gene.start-1, gene.end, gene.strand, name=gene.ensembl_id, original_base=1)
                        all_regions.append(r)

                elif args.exons:
                    exon_ids = {}
                    transcripts = gtf_db.get_transcripts_simple(args.database)
                    for i, transcript in enumerate(transcripts):
                        for ensembl_id, exon in transcript.exons.iteritems():
                            if ensembl_id not in exon_ids:
                                r = g2g.Region(exon.seqid, exon.start-1, exon.end, exon.strand, name=exon.ensembl_id, original_base=1)
                                all_regions.append(r)
                                exon_ids[ensembl_id]=1

            fasta.extract(args.fasta, all_regions, output=None, reverse=reverse, complement=complement, raw=args.raw)

    except KeyboardInterrupt as ki:
        g2g.exit("Handled Keyboard Interrupt")
    except exceptions.G2GValueError as ve:
        g2g.exit(ve.msg)
    except exceptions.G2GBedError as be:
        g2g.exit(be.msg)
    except exceptions.G2GRegionError as le:
        g2g.exit(le.msg)
    except exceptions.G2GFastaError as fe:
        g2g.exit(fe.msg)


def command_fasta_patch(raw_args, prog=None):
    """
    Convert a Fasta file into a patched (SNP'd) Fasta file.

    Usage: patch [-options] -i <Fasta file> -v <VCI file>

    Required Parameters:
        -i, --input <Fasta file>         Fasta file to patch
        -v, --vci <VCI file>             input VCI file

    Optional Parameters:
        -b, --bed <BED file>             BED file (cannot be used with -l)
        --bgzip                          compress and index output (can only be used with -o)
        -o, --output <Output file>       Fasta file to output to
        -n, --num_processes <number>     the number of processes to use, defaults to the number of cores
        -r, --region <region>            seqid:start-end (cannot be used with -b)
        --reverse                        reverse the VCI file

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    Note: --bgzip can be potentially slow depending on the size of the file

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_fasta_patch.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")
    parser.add_argument("-v", "--vci", dest="vci", metavar="VCI_File")

    # optional
    parser.add_argument("--bgzip", dest="bgzip", action='store_true')
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-n", "--numprocesses", type=int, dest="numprocesses", metavar="number_of_processes")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")
    parser.add_argument("-r", "--region", dest="region", metavar="chr1:100000-200000")
    parser.add_argument("--reverse", dest="reverse", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    if not args.fasta:
        g2g.exit("No input Fasta file was specified.", parser)

    if not args.vci:
        g2g.exit("No VCI file was specified.", parser)

    if not args.output and args.bgzip:
        g2g.exit("--bgzip can only be used with the -o option", parser)

    all_locations = None

    if args.region and args.bed_file:
        g2g.app_exit("Please use either a location or a BED file.", parser)
    elif args.region and not args.bed_file:
        all_locations = [g2g.parse_region(args.region)]
    elif not args.region and args.bed_file:
        bed_file = bed.BED(args.bed_file)
        all_locations = []
        for bed_rec in bed_file:
            if bed_file.current_line_is_bed:
                strand = bed_rec.strand if bed_rec.strand else '+'
                all_locations.append(g2g.Region(bed_rec.chrom, bed_rec.start, bed_rec.end, strand, name=bed_rec.name))

    try:
        fasta_patch.process(args.fasta, args.vci, all_locations, args.output, args.bgzip, args.reverse, args.numprocesses)
    except exceptions.KeyboardInterruptError as ki:
        g2g.exit(ki, parser)
    except exceptions.G2GValueError as e:
        g2g.exit(e, parser)
    except exceptions.G2GVCFError as e:
        g2g.exit(e, parser)
    except exceptions.G2GError as e:
        g2g.exit(e)


def command_fasta_transform(raw_args, prog=None):
    """
    Extract a sequence from a Fasta file give a location, and replace the from
    newly formatted G2G file.

    Usage: transform [-options] [-options] -i <Fasta file> -v <VCI file>

    Required Parameters:
        -i, --input <Fasta file>         Fasta file to extract sequence from
        -v, --vci <VCI file>             input VCI file

    Optional Parameters:
        -b, --bed <BED file>             BED file (cannot be used with -l)
        --bgzip                          compress and index output (can only be used with -o)
        -l, --location <location>        seqid:start-end (cannot be used with -b)
        -o, --output <Output file>       Fasta file to output to
        -n, --num_processes <number>     the number of processes to use, defaults to the number of cores
        -r, --reverse                    reverse the chain file

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    Note: --bgzip can be potentially slow depending on the size of the file

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_fasta_transform.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")
    parser.add_argument("-v", "--vci", dest="vci", metavar="VCI_File")

    # optional
    parser.add_argument("--bgzip", dest="bgzip", action='store_true')
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-n", "--numprocesses", type=int, dest="numprocesses", metavar="number_of_processes")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")
    parser.add_argument("-r", "--region", dest="region", metavar="chr1:100000-200000")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    if not args.fasta:
        g2g.exit("No input Fasta file was specified.", parser)

    if not args.vci:
        g2g.exit("No VCI file was specified.", parser)

    if not args.output and args.bgzip:
        g2g.exit("--bgzip can only be used with the -o option", parser)

    try:
        all_locations = None

        if args.region and args.bed_file:
            g2g.app_exit("Please use either a location or a BED file.", parser)
        elif args.region and not args.bed_file:
            all_locations = [g2g.parse_region(args.region)]
        elif not args.region and args.bed_file:
            bed_file = bed.BED(args.bed_file)
            all_locations = []
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else '+'
                    all_locations.append(g2g.Region(bed_rec.chrom, bed_rec.start, bed_rec.end, strand))

        fasta_transform.process(args.fasta, args.vci, all_locations, args.output, args.bgzip, False, args.numprocesses, False)
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GChainFileError as e:
        g2g.exit("Chain File error")
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GRegionError as e:
        g2g.exit("Region error")
    except exceptions.G2GFastaError as e:
        g2g.exit("Fasta error")
    except exceptions.G2GError as e:
        g2g.exit(e)


def command_parse_region(raw_args, prog=None):
    """
    Parse a region from the command line.  Useful for verifying regions.

    Usage: parse [-options] -l <location>

    Required Parameters:
        -r, --region <region>            seqid:start-end

    Optional Parameters:
        -b, --base <0|1>                 base (0 or 1)
        -n, --name <name>                name of region

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_parse_region.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-r", "--region", dest="region", metavar="chr1:100000-200000")

    # optional
    parser.add_argument("-b", "--base", dest="base", metavar="0", default=None)
    parser.add_argument("-n", "--name", dest="name", default=None)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    try:
        if args.region:
            LOG.debug("Input: {0}".format(args.region))

            if args.name:
                LOG.debug("-n '{}'".format(args.name))

            base = 1

            if args.base:
                LOG.debug("-b {}".format(args.base))
                base = args.base

            region = g2g.parse_region(args.region, base=base, name=args.name)
            LOG.info(region)
            LOG.info("Seq ID: {}".format(region.seq_id))
            LOG.info("Start: {}".format(region.start))
            LOG.info("_Start: {}".format(region._start))
            LOG.info("End: {}".format(region.end))
            LOG.info("Name: {}".format(region.name))
            LOG.info("Base: {}".format(region.original_base))
            LOG.info("Display Start: {}".format(region.get_start()))

    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as ve:
        g2g.exit(ve.msg)
    except exceptions.G2GBedError as be:
        g2g.exit(be.msg)
    except exceptions.G2GRegionError as le:
        g2g.exit(le.msg)
    except exceptions.G2GFastaError as fe:
        g2g.exit(fe.msg)


def command_fastaformat(raw_args, prog=None):
    """
    Format Fasta file

    Usage: fastaformat [-options] -f <Fasta file>

    Required Parameters:
        -f, --fasta <fasta_file>         Fasta file to format

    Optional Parameters:
        -l, --length <line_length>       the length of the line, defaults to 60
        -o, --output <output_file>       Fasta file to create
        -s, --seqids <seqids>            comma seperated list of seqids

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_fastaformat.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-f", "--fasta", dest="fasta", metavar="Fasta_File")

    # optional
    parser.add_argument("-l", "--length", dest="line_length", metavar="line_length", default=60, type=int)
    parser.add_argument("-o", "--output", dest="output", metavar="fasta_file")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if len(raw_args) == 0 or args.help:
        g2g.exit("", parser)

    if not args.fasta:
        g2g.exit("No Fasta file was specified.", parser)

    try:
        fasta.reformat(args.fasta, args.output, args.line_length)
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as e:
        g2g.exit(e, parser)
    except exceptions.G2GVCFError as e:
        g2g.exit(e, parser)


def command_gtf2db(raw_args, prog=None):
    """
    Convert a GTF file to a G2G DB file

    Usage: gtf2db [-options] -i <GTF file> -o <G2G DB file>

    Required Parameters:
    -i, --input <GTF file>           GTF file to parse
    -o, --output <G2G DB file>       G2G DB file to create

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_gtf2db.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="GTF_file")
    parser.add_argument("-o", "--output", dest="output", metavar="DB_file")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if args.help:
        g2g.exit("", parser)

    if not args.input:
        g2g.exit("No GTF file was specified.", parser)

    if not args.output:
        g2g.exit("No output GTG DB file was specified.", parser)

    try:
        gtf_db.gtf2db(args.input, args.output)
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as e:
        g2g.exit(e, parser)
    except exceptions.G2GError as e:
        g2g.exit(e, parser)


def command_vciquery(raw_args, prog=None):
    """
    Query a VCI file

    Usage: vciquery [-options] -i <VCI file>

    Required Parameters:
        -i, --input <VCI file>           VCI file to query
        -r, --region                     Region to query

    Optional Parameters:
        -f, --fasta <Fasta file>         Fasta file

    Help Parameters:
        -h, --help                       print the help and exit
        -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_vciquery.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-v", "--vci", dest="vci", metavar="VCI_file")
    parser.add_argument("-r", "--region", dest="region", metavar="chr1:100000-200000")

    # optional
    parser.add_argument("-f", "--fasta", dest="fasta", metavar="Fasta_file")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    g2g.configure_logging(args.debug)

    if args.help:
        g2g.exit("", parser)

    if not args.vci:
        g2g.exit("No VCI file was specified.", parser)

    if not args.region:
        g2g.exit("No region was specified.", parser)

    try:
        region = g2g.parse_region(args.region, 1)
        vci.vci_query(args.vci, region, args.fasta)
    except KeyboardInterrupt as ki:
        LOG.debug(ki)
    except exceptions.G2GValueError as e:
        g2g.exit(e, parser)
    except exceptions.G2GError as e:
        g2g.exit(e, parser)

# -*- coding: utf-8 -*-

import argparse
import sys

from .bed import BED
from .exceptions import G2GError, G2GBAMError, G2GBedError, G2GChainFileError, G2GVCFError, G2GValueError, G2GLocationError, G2GFastaError
from .g2g_utils import Location, app_exit, configure_logging, get_logger, parse_location
from .g2g import file_convert, gtf2chain, db2chain, offset2chain, fasta_extract_exons, fasta_extract_transcripts, dbfetch, fasta_transform
from .gtf_db import gtf2db, get_genes_simple
from .fasta import extract as fasta_extract
from .fasta import extract_id as fasta_extract_id
from .fasta_patch import fasta_patch
from .vcf2chain import vcf2chain
LOG = get_logger()


def command_convert(raw_args, prog=None):
    """
    Convert a file

    Usage: convert [-options] -c <Chain file> -i <input file>

    Required Parameters:
    -c, --chain <Chain file>         input Chain file used in conversion
    -i, --input <input file>         input file to convert
    -o, --output <output file>       file to create

    Optional Parameters:
    -f, --format <bam|sam|gtf|bed>   file format
    -r, --reverse                    reverse the Chain file

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
        print message
        print command_convert.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-c", "--chain", dest="chain", metavar="Chain_File")
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")

    # optional
    parser.add_argument("-f", "--format", dest="format", metavar='bam|sam|gtf|bed', choices=['bam', 'sam', 'gtf', 'bed'])
    parser.add_argument("-r", "--reverse", dest="reverse", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.chain:
        app_exit("No Chain file was specified.", parser)

    if not args.input:
        app_exit("No input file was specified.", parser)

    if not args.output:
        app_exit("No output file was specified.", parser)

    try:

        file_convert(args.chain, args.input, args.output, args.format, args.reverse)

    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, ve:
        app_exit(ve, parser)
    except G2GChainFileError, cfe:
        app_exit(cfe, parser)
    except G2GBAMError, e:
        app_exit(e, parser)


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
        print message
        print command_gtf2db.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="GTF_file")
    parser.add_argument("-o", "--output", dest="output", metavar="DB_file")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.input:
        app_exit("No GTF file was specified.", parser)

    if not args.output:
        app_exit("No output GTG DB file was specified.", parser)

    try:
        gtf2db(args.input, args.output)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GChainFileError, e:
        app_exit(e, parser)
    except G2GError, e:
        app_exit(e, parser)


def command_gtf2chain(raw_args, prog=None):
    """
    Convert a GTF file to a new chain file.

    Usage: gtf2chain [-options] -c <Chain file> -i <GTF file>

    Required Parameters:
    -c, --chain <Chain file>         input Chain file used in conversion
    -i, --input <GTF file>           GTF file to parse
    -o, --output <Chain file>        Chain file to create

    Optional Parameters:
    -g, --genes                      Create Gene Chain File, defaults to Transcript chain file

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """
    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_gtf2chain.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-c", "--chain", dest="chain", metavar="Chain_File")
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Chain_File")

    # optional
    parser.add_argument("-g", "--genes", dest="genes", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.chain:
        app_exit("No from Chain file was specified.", parser)

    if not args.input:
        app_exit("No GTF file was specified.", parser)

    if not args.output:
        app_exit("No output Chain file was specified.", parser)

    try:
        gtf2chain(args.chain, args.input, args.output, args.genes)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GChainFileError, e:
        app_exit(e, parser)
    except G2GError, e:
        app_exit(e, parser)


def command_db2chain(raw_args, prog=None):
    """
    Convert a G2G DB file to a Chain file

    Usage: db2chain [-options] -c <Chain file> -i <G2G DB file>

    Required Parameters:
    -c, --chain <Chain file>         input Chain file used in conversion
    -i, --input <G2G DB file>        G2F DB file
    -o, --output <Chain file>        Chain file to create

    Optional Parameters:
    -g, --genes                      Create Gene Chain File, defaults to Transcript chain file

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_db2chain.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-c", "--chain", dest="chain", metavar="Chain_File")
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Chain_File")

    # optional
    parser.add_argument("-g", "--genes", dest="genes", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.chain:
        app_exit("No input Chain file was specified.", parser)

    if not args.input:
        app_exit("No G2G DB file was specified.", parser)

    if not args.output:
        app_exit("No output Chain was specified.", parser)

    try:
        db2chain(args.chain, args.input, args.output, args.genes)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GChainFileError, e:
        app_exit(e, parser)


def command_vcf2chain(raw_args, prog=None):
    """
    Convert a VCF file to a Chain file

    Usage: vcf2chain [-options] -i <VCF file> -f <Fasta file> -s <strain>

    Required Parameters:
    -f, --fasta <Fasta file>         Fasta file to parse, BGZIP file supported
    -i, --input <input file>         VCF file to convert
    -o, --output <output file>       Chain file to create
    -s, --strain <strain>            name of strain column in VCF file

    Optional Parameters:
    --diploid                        ignore hetorzygotes, create 2 files
    --keep                           keep track of VCF lines that cannot be converted to Chain file
    --pass                           use only VCF lines that have a PASS for the filter value
    --quality                        filter on quality, FI=PASS

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_vcf2chain.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-f", "--fasta", dest="fasta", metavar="Fasta_File")
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Chain_File")
    parser.add_argument("-s", "--strain", dest="strain", metavar="strain")

    # optional
    parser.add_argument("--diploid", dest="diploid", action='store_true')
    parser.add_argument("--keep", dest="keep", action='store_true')
    parser.add_argument("--pass", dest="passed", action='store_true')
    parser.add_argument("--quality", dest="quality", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.input:
        app_exit("No VCF file was specified.", parser)

    if not args.fasta:
        app_exit("No Fasta file was specified.", parser)

    if not args.output:
        app_exit("No Chain file was specified.", parser)

    if not args.strain:
        app_exit("No strain was specified.", parser)

    try:
        vcf2chain(args.input, args.fasta, args.strain, args.output, args.keep, args.passed, args.quality, args.diploid)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GVCFError, e:
        app_exit(e, parser)


def command_offset2chain(raw_args, prog=None):
    """
    Convert Seqnature offset files to a Chain file

    Usage: offset2chain [-options] -f <length_file_of_from-genome> -t <length_file_of_to-genome>

    Required Parameters:
      -f, --from <length_file_of_from-genome>  Chromosome length file of From-genome
      -o, --output <output file>  Chain file to create
      -t, --to   <length_file_of_to-genome>    Chromosome length file of To-genome

    Optional Parameters:
      None

    Help Parameters:
      -h, --help                  Print the help and exit
      -d, --debug                 Turn debugging on, list multiple times for more messages

    Note:
      Consider using the following awk command if you don't have chromosome length files.
      $ awk '/^>/ { if(seqlen) { print seqlen }; { printf "%s\t", substr($1,2) }; seqlen=0; next; } \
        { seqlen=seqlen+length($0) } END { print seqlen }' ${genome_fasta} > ${outfile}

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_offset2chain.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-f", "--from", dest="fromchrom", metavar="Chromosome_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Chain_File")
    parser.add_argument("-t", "--to", dest="tochrom", metavar="Chromosome_File")

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.fromchrom:
        app_exit("No from file was specified.", parser)

    if not args.tochrom:
        app_exit("No from to was specified.", parser)

    if not args.output:
        app_exit("No output file was specified.", parser)

    try:
        offset2chain(args.fromchrom, args.tochrom, args.output)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GChainFileError, ce:
        app_exit(ce, parser)


def command_fasta_extract(raw_args, prog=None):
    """
    Extract a sequence from a Fasta file give a location.

    Usage: extract [-options] -i <Fasta file> [-l <location> | -b <BED file> | -db <Database> | -id <id>]

    Required Parameters:
    -i, --input <Fasta file>         Fasta file to extract sequence from, BGZIP files supported
    -l, --location <location>        seqid:start:end
        - OR -
    -b, --bed <BED file>             BED file
        - OR -
    -db, --database <DB file>        Database file
        - OR -
    -id, --identifier <identifier>   Fasta identifier

    Optional Parameters:
    -c, --complement                 return complement of sequence
    -e, --exons                      For use with -db, defaults to transcript extraction
    -g, --genes                      For use with -db, defaults to transcript extraction
    --raw                            show just sequence
    -r, --reverse                    return reverse of sequence
    -rc, --reversecomplement         return reverse complement of sequence
    -t, --transcripts                For use with -db, defaults to transcript extraction

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
        print message
        print command_fasta_extract.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-id", "--identifier", dest="id", metavar="identifier")
    parser.add_argument("-l", "--location", dest="location", metavar="chr1:100000-200000")
    parser.add_argument("-db", "--database", dest="database", metavar="Database")

    # optional
    parser.add_argument("-c", "--complement", dest="complement", action="store_true", default=False)
    parser.add_argument("-e", "--exons", dest="exons", action="store_true", default=False)
    parser.add_argument("-g", "--genes", dest="genes", action="store_true", default=False)
    parser.add_argument("-r", "--reverse", dest="reverse", action="store_true", default=False)
    parser.add_argument("--raw", dest="raw", action="store_true", default=False)
    parser.add_argument("-rc", "--reversecomplement", dest="reversecomplement", action="store_true", default=False)
    parser.add_argument("-t", "--transcripts", dest="transcripts", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.fasta:
        app_exit("No Fasta file was specified.", parser)

    flag_count = 0
    flag_count += (1 if args.complement else 0)
    flag_count += (1 if args.reverse else 0)
    flag_count += (1 if args.reversecomplement else 0)

    if flag_count > 1:
        app_exit("Please specify only one of: -c, -r, -rc", parser)

    reverse = (args.reverse or args.reversecomplement)
    complement = (args.complement or args.reversecomplement)

    if args.database:
        # allow only either exon, gene, or transcript extraction, bot all
        output_count = 0
        output_count += (1 if args.genes else 0)
        output_count += (1 if args.transcripts else 0)
        output_count += (1 if args.exons else 0)

        if output_count == 0:
            app_exit("Please specify one of: -e, -g, or -t", parser)
        elif output_count > 1:
            app_exit("Please specify only one of: -e, -g, or -t", parser)

    location_count = 0
    location_count += (1 if args.location else 0)
    location_count += (1 if args.bed_file else 0)
    location_count += (1 if args.database else 0)
    location_count += (1 if args.id else 0)

    if location_count > 1:
        app_exit("Please specify only one of: -l, -b, -id, or -db", parser)

    all_locations = []

    try:
        if args.location:
            all_locations.append(parse_location(args.location, 1))

        elif args.bed_file:
            if flag_count == 1:
                app_exit("Options -c, -r, -rc cannot be used with BED file.", parser)

            bed_file = BED(args.bed_file)
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else '+'
                    all_locations.append(Location(bed_rec.chrom, bed_rec.start, bed_rec.end, strand, bed_rec.name, 0))

        elif args.database:
            if flag_count == 1:
                app_exit("Options -c, -r, -rc cannot be used with Database.", parser)

            if args.genes:
                LOG.info("Extracting genes from database...")
                genes = get_genes_simple(args.database)

                LOG.debug("Genes extracted, transforming to locations")
                for gene in genes:
                    l = Location(gene.seqid, gene.start-1, gene.end, gene.strand, name=gene.ensembl_id, original_base=1)
                    all_locations.append(l)

                LOG.info("Extracting genes sequences...")
                fasta_extract(args.fasta, all_locations, None, reverse=reverse, complement=complement, raw=args.raw)
                LOG.info("Extracting genes done")

            elif args.transcripts:
                fasta_extract_transcripts(args.fasta, args.database, None, raw=args.raw)

            elif args.exons:
                fasta_extract_exons(args.fasta, args.database, None, raw=args.raw)
            else:
                LOG.error("Nothing to do")

            return

        elif args.id:
            fasta_extract_id(args.fasta, args.id, output=None, reverse=reverse, complement=complement, raw=args.raw)
            return

        fasta_extract(args.fasta, all_locations, output=None, reverse=reverse, complement=complement, raw=args.raw)
    except KeyboardInterrupt, ki:
        LOG.info("Handled Keyboard Interrupt")
    except G2GValueError, ve:
        app_exit(ve.msg)
    except G2GBedError, be:
        app_exit(be.msg)
    except G2GLocationError, le:
        app_exit(le.msg)
    except G2GFastaError, fe:
        app_exit(fe.msg)


def command_fasta_transform(raw_args, prog=None):
    """
    Extract a sequence from a Fasta file give a location, and replace the from
    newly formatted chain file.

    Usage: transform [-options] -c <Chain file> -i <Fasta file> [-l <location> | -b <BED file>]

    Required Parameters:
    -c, --chain <Chain file>         input Chain file
    -i, --input <Fasta file>         Fasta file to extract sequence from
    -l, --location <location>        seqid:start:end
        - OR -
    -b, --bed <BED file>             BED file
    -o, --output <Output file>       file to output to

    Optional Parameters:
    --bgzip                          compress and index output
    -c, --complement                 return complement of sequence
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
        print message
        print command_fasta_transform.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-c", "--chain", dest="chain", metavar="Chain_File")
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-l", "--location", dest="location", metavar="chr1:100000-200000")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")

    # optional
    parser.add_argument("--bgzip", dest="bgzip", action='store_true')
    parser.add_argument("--diploid", dest="diploid", action='store_true')
    parser.add_argument("-r", "--reverse", dest="reverse", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.chain:
        app_exit("No Chain file was specified.", parser)

    if not args.fasta:
        app_exit("No input Fasta file was specified.", parser)

    if not args.output:
        app_exit("No output Fasta file was specified.", parser)

    try:
        all_locations = None

        if args.location and args.bed_file:
            app_exit("Please use either a location or a BED file.", parser)
        elif args.location and not args.bed_file:
            all_locations = [parse_location(args.location)]
        elif not args.location and args.bed_file:
            bed_file = BED(args.bed_file)
            all_locations = []
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else '+'
                    all_locations.append(Location(bed_rec.chrom, bed_rec.start, bed_rec.end, strand))

        fasta_transform(args.fasta, args.chain, all_locations, args.output, args.bgzip, args.reverse)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GChainFileError, le:
        app_exit("Chain File error")
    except G2GValueError, le:
        app_exit(le.msg)
    except G2GLocationError, le:
        app_exit("Location error")
    except G2GFastaError, fe:
        app_exit("Fasta error")


def command_fasta_patch(raw_args, prog=None):
    """
    Convert a Fasta file into a patched (SNP'd) Fasta file

    Usage: patch [-options] -i <Fasta file> -v <VCF file> -s <strain> -o <Fasta File>

    Required Parameters:
    -i, --input <Fasta file>         Fasta file to patch
    -o, --output <Output file>       Fasta file to output to
    -s, --strain <strain>            name of strain column in VCF file
    -v, --vcf <input file>           VCF file to parse SNPs

    Optional Parameters:
    --bgzip                          compress and index output
    --diploid                        ignore hetorzygotes, create 2 files
    -n, --num_processes <number>     the number of processes to use, defaults to the number of cores
    -p, --pass                       use only VCF lines that have a PASS for the filter value
    --quality                        filter on quality, FI=PASS

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
        print message
        print command_fasta_patch.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="fasta", metavar="FASTA_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")
    parser.add_argument("-s", "--strain", dest="strain", metavar="strain")
    parser.add_argument("-v", "--vcf", dest="vcf", metavar="Input_File")

    # optional
    parser.add_argument("--bgzip", dest="bgzip", action='store_true')
    parser.add_argument("--diploid", dest="diploid", action='store_true')
    parser.add_argument("-n", "--numprocesses", type=int, dest="numprocesses", metavar="number_of_processes")
    parser.add_argument("-p", "--pass", dest="passonly", action='store_true')
    parser.add_argument("--quality", dest="quality", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.fasta:
        app_exit("No input Fasta file was specified.", parser)

    if not args.vcf:
        app_exit("No VCF file was specified.", parser)

    if not args.strain:
        app_exit("No strain was specified.", parser)

    if not args.output:
        app_exit("No output Fasta file was specified.", parser)

    try:
        fasta_patch(args.fasta, args.vcf, args.strain, args.output, args.bgzip, args.numprocesses, args.passonly, args.quality, args.diploid)
    except KeyboardInterrupt, ki:
        app_exit(ki, parser)
    except G2GValueError, e:
        app_exit(e, parser)
    except G2GVCFError, e:
        app_exit(e, parser)
    except G2GError, e:
        sys.exit(1)


def command_dbfetch(raw_args, prog=None):
    """
    Retrieve information from a G2G Database

    Usage: dbfetch [-options] -i <G2G DB file> [-l <location> | -b <BED file>]

    Required Parameters:
    -i, --input <G2G DB file>        G2F DB file
    -l, --location <location>        seqid:start:end
        - OR -
    -b, --bed <BED file>             BED file
        - OR -
    --id <Ensembl ID>                Ensembl ID


    Optional Parameters:
    -c, --inclusive                  inclusive or exclusive
    -e, --exons                      output exon information
    -g, --genes                      output gene information
    -o, --output <Output file>       redirect output to a file
    -t, --transcripts                output transcript information

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_dbfetch.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-b", "--bed", dest="bed_file", metavar="BED_File")
    parser.add_argument("-l", "--location", dest="location", metavar="chr1:100000-200000")
    parser.add_argument("--id", dest="id", metavar="ENSMUSG00000000000")

    # optional
    parser.add_argument("-c", "--inclusive", dest="inclusive", action="store_true", default=False)
    parser.add_argument("-e", "--exons", dest="exons", action="store_true", default=False)
    parser.add_argument("-g", "--genes", dest="genes", action="store_true", default=False)
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")
    parser.add_argument("-t", "--transcripts", dest="transcripts", action="store_true", default=False)

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    if not args.input:
        app_exit("No G2G DB file was specified.", parser)

    location_count = 0
    location_count += (1 if args.location else 0)
    location_count += (1 if args.bed_file else 0)
    location_count += (1 if args.id else 0)

    if location_count != 1:
        app_exit("Please specify either -l or -b or --id", parser)

    if (args.genes and args.transcripts and args.exons) or (not args.genes and not args.transcripts and not args.exons):
        args.genes = True
        args.transcripts = True
        args.exons = True

    all_locations = []

    try:
        if args.id:
            dbfetch(args.input, args.output, ids=[args.id], display_genes=args.genes, display_transcripts=args.transcripts, display_exons=args.exons)
        else:
            if args.location:
                all_locations.append(parse_location(args.location, 1))

            elif args.bed_file:
                bed_file = BED(args.bed_file)
                for bed_rec in bed_file:
                    if bed_file.current_line_is_bed:
                        strand = bed_rec.strand if bed_rec.strand else '+'
                        all_locations.append(Location(bed_rec.chrom, bed_rec.start, bed_rec.end, strand, 0))

            dbfetch(args.input, args.output, locations=all_locations, display_genes=args.genes, display_transcripts=args.transcripts, display_exons=args.exons)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, ve:
        app_exit(ve.msg)
    except G2GBedError, be:
        app_exit(be.msg)
    except G2GLocationError, le:
        app_exit(le.msg)
    except G2GFastaError, fe:
        app_exit(fe.msg)


def command_parse_location(raw_args, prog=None):
    """
    Parse a location from the command line.  Useful for verifying locations.

    Usage: parse [-options] -l <location>

    Required Parameters:
    -l, --location <location>        seqid:start:end

    Optional Parameters:
    None

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message):
        print message
        print command_fasta_extract.__doc__
        sys.exit()

    parser.error = print_message

    # required
    parser.add_argument("-l", "--location", dest="location", metavar="chr1:100000-200000")

    # optional

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    configure_logging(args.debug)

    if args.help:
        app_exit("", parser)

    try:
        if args.location:
            LOG.debug("Input: {0}".format(args.location))
            loc = parse_location(args.location, 1)
            LOG.debug(loc.__repr__())
            print str(loc)

    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except G2GValueError, ve:
        app_exit(ve.msg)
    except G2GBedError, be:
        app_exit(be.msg)
    except G2GLocationError, le:
        app_exit(le.msg)
    except G2GFastaError, fe:
        app_exit(fe.msg)

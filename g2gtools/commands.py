# standard library imports
# none

# 3rd party library imports
import click

# local library imports
from g2gtools import __version__ as version
import g2gtools.bcsam as bcsam
import g2gtools.bed as g2g_bed
import g2gtools.exceptions as exceptions
import g2gtools.fasta as g2g_fasta
import g2gtools.fasta_patch as g2g_fasta_patch
import g2gtools.fasta_transform as g2g_fasta_transform
import g2gtools.g2g as g2g
import g2gtools.g2g_utils as g2g_utils
import g2gtools.gtf as gtf
import g2gtools.gtf_db as gtf_db
import g2gtools.vci as g2g_vci
import g2gtools.vcf2vci as g2g_vcf2vci

logo_text = r"""

        ___       _              _
       |__ \     | |            | |
   __ _   ) |__ _| |_ ___   ___ | |___
  / _` | / // _` | __/ _ \ / _ \| / __|
 | (_| |/ /| (_| | || (_) | (_) | \__ \\
  \__, |____\__, |\__\___/ \___/|_|___/   v""" + version + """
   __/ |     __/ |
  |___/     |___/

"""


@click.group()
@click.pass_context
def cli(ctx):
    # ensure that ctx.obj exists and is a dict
    ctx.ensure_object(dict)


# #############################################################################
#
# convert
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--input-file",
    help="Input file to convert to new coordinates",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-c",
    "--vci",
    help="VCI file",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-f",
    "--file-format",
    help="Input file format",
    required=True,
    type=click.Choice(["BAM", "SAM", "GTF", "BED"], case_sensitive=False)
)
@click.option(
    "-o",
    "--output-file",
    help="Name of output file",
    type=click.Path(
        exists=False, dir_okay=False, writable=True, resolve_path=True
    )
)
@click.option(
    "-r",
    "--reverse",
    default=False,
    help="Reverse the direction of the conversion",
    is_flag=True
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def convert(ctx, input_file, vci, file_format, output_file, reverse, verbose):
    """
    Convert coordinates of BAM|SAM|GTF|GFF|BED files
    """
    try:
        if file_format in ["BAM", "SAM"]:
            bcsam.convert_bcsam_file(
                vci_file=vci,
                bcsam_file_name_in=input_file,
                bcsam_file_name_out=output_file,
                reverse=reverse,
                debug_level=verbose
            )
        elif file_format in ["GTF"]:
            gtf.convert_gtf_file(
                vci_file=vci,
                gtf_file_name_in=input_file,
                gtf_file_name_out=output_file,
                reverse=reverse,
                debug_level=verbose
            )
        elif file_format in ["BED"]:
            g2g_bed.convert_bed_file(
                vci_file=vci,
                bed_file_name_in=input_file,
                bed_file_name_out=output_file,
                reverse=reverse,
                debug_level=verbose,
            )
        elif file_format in ["GFF"]:
            gtf.convert_gff_file(
                vci_file=vci,
                gff_file_name_in=input_file,
                gff_file_name_out=output_file,
                reverse=reverse,
                debug_level=verbose,
            )
        else:
            raise exceptions.G2GValueError(f"Unknown file format: {format}")
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GChainFileError as e:
        g2g.exit(e.msg)
    except exceptions.G2GBAMError as e:
        g2g.exit(e.msg)
    except exceptions.G2GBedError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# vcf2vci
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--vcf",
    multiple=True,
    help="VCF file name",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-f",
    "--fasta",
    help="Fasta file matching VCF information",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-s",
    "--strain",
    help="Name of strain/sample (column in VCF file)",
    type=str,
    required=True
)
@click.option(
    "-o",
    "--output",
    help="VCI file name to create",
    required=True,
    type=click.Path(
        exists=False, dir_okay=False, writable=True, resolve_path=True
    )
)
@click.option(
    "-p",
    "--num-processes",
    help="The number of processes to use, defaults to the number of cores",
    type=int
)
@click.option(
    "--diploid",
    default=False,
    help="Create diploid VCI file",
    is_flag=True
)
@click.option(
    "--keep",
    default=False,
    help="Keep track of VCF lines that could not be converted to VCI file",
    is_flag=True
)
@click.option(
    "--pass",
    "passed",
    default=False,
    help="Use only VCF lines that have a PASS for the filter value",
    is_flag=True
)
@click.option(
    "--quality",
    default=False,
    help="Filter on quality, FI=PASS",
    is_flag=True
)
@click.option(
    "--no-bgzip",
    default=False,
    help="DO NOT compress and index output",
    is_flag=True
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def vcf2vci(ctx, vcf, fasta, strain, output, keep, passed, quality, diploid,
            num_processes, no_bgzip, verbose):
    """
    Create VCI file from VCF file(s).
    """
    try:
        g2g_vcf2vci.process(
            vcf_files=vcf,
            fasta_file=fasta,
            output_file=output,
            strain=strain,
            vcf_keep=keep,
            passed=passed,
            quality=quality,
            diploid=diploid,
            num_processes=num_processes,
            bgzip=(not no_bgzip),
            debug_level=verbose
        )
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GVCFError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# extract
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--fasta",
    help="Fasta file to extract from",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-b",
    "--bed",
    help="BED file name",
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-db",
    "--database",
    help="Database file name, use with --genes, --transcripts, --exons",
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "--genes",
    default=False,
    help="Extract genes from --database",
    is_flag=True
)
@click.option(
    "--transcripts",
    default=False,
    help="Extract transcripts from --database",
    is_flag=True
)
@click.option(
    "--exons",
    default=False,
    help="Extract exons from --database",
    is_flag=True
)
@click.option(
    "-id",
    "--identifier",
    help="Fasta identifier",
    type=str
)
@click.option(
    "-r",
    "--region",
    help="Region to extract in chromosome:start-end format",
    type=str
)
@click.option(
    "-c",
    "--vci",
    help="VCI file to use",
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-R",
    "--vci-reverse",
    default=False,
    help="Reverse the direction of the VCI file",
    is_flag=True
)
@click.option(
    "--complement",
    default=False,
    help="Complement the extracted sequence",
    is_flag=True
)
@click.option(
    "--reverse",
    default=False,
    help="Reverse the extracted sequence",
    is_flag=True
)
@click.option(
    "--reverse-complement",
    default=False,
    help="Reverse-complement the extracted sequence",
    is_flag=True
)
@click.option(
    "--raw",
    default=False,
    help="Shows just the extracted sequences",
    is_flag=True
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def extract(ctx, fasta, bed, database, genes, transcripts, exons, identifier,
            region, vci, vci_reverse, complement, reverse, reverse_complement,
            raw, verbose):
    """
    Extract subsequence from a fasta file given a region

    Note:
        Locations specified on command line are 1-based coordinates.
        Locations specified via BED file are 0-based coordinates.
    """
    flag_count = 0
    flag_count += 1 if complement else 0
    flag_count += 1 if reverse else 0
    flag_count += 1 if reverse_complement else 0

    if flag_count > 1:
        g2g.exit(
            "Please specify only one of: --complement, --reverse, "
            "--reverse-complement",
        )

    reverse = reverse or reverse_complement
    complement = complement or reverse_complement

    if database:
        # allow only either exon, gene, or transcript extraction, bot all
        output_count = 0
        output_count += 1 if genes else 0
        output_count += 1 if transcripts else 0
        output_count += 1 if exons else 0

        if output_count != 1:
            g2g.exit(
                "Please specify one of: --exons, --genes, or --transcripts"
            )

    location_count = 0
    location_count += 1 if region else 0
    location_count += 1 if bed else 0
    location_count += 1 if database else 0
    location_count += 1 if identifier else 0

    if location_count > 1:
        g2g.exit("Please specify only one of: -r, -b, -id, or -db")

    logger = g2g.get_logger(verbose)

    # parse the regions we are extracting
    all_regions = []

    try:
        fasta_file = g2g_utils.check_file(fasta)
        fasta_file = g2g_fasta.FastaFile(fasta_file)

        # get the regions if there are some

        if region:
            # simple region
            logger.debug(f"ARGS.REGION={region}")
            region = g2g.parse_region(region, base=1)

            logger.debug(f"--> start = {region.start}")
            logger.debug(f"--> _start = {region._start}")
            logger.debug(f"--> get_start() = {region.get_start()}")

            logger.debug(f"REGION PARSED={region}")
            all_regions.append(region)

            # TEMPORARY OVERRIDE
            # all_regions = [g2g.parse_region(args.region)]

        elif bed:
            # bed file
            bed_file = bed.BED(bed)
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else "+"
                    logger.debug(bed_rec)
                    all_regions.append(
                        g2g.Region(
                            bed_rec.chrom,
                            bed_rec.start,
                            bed_rec.end,
                            strand,
                            bed_rec.name,
                            0,
                        )
                    )

        # on the fly
        if vci:
            # let's eliminate what isn't supported

            if fasta_file.is_diploid():
                g2g.exit("Diploid Fasta files are not supported with -v, --vci")
            elif identifier:
                g2g.exit(
                    "Options -v, --vci cannot be used with --id, --identifier"
                )
            elif database:
                g2g.exit(
                    "Options -v, --vci cannot be used with --db, --database"
                )
            elif flag_count >= 1:
                g2g.exit(
                    "Options --complement, --reverse, --reversecomplement "
                    "cannot be used with VCI."
                )

            if region:
                # overridden here, fasta.extract expects 1 base
                all_regions = [g2g.parse_region(region)]

            g2g_fasta_transform.process(
                fasta_file_name_in=fasta,
                vci_file_name=vci,
                regions=all_regions,
                fasta_file_name_out=None,
                bgzip=False,
                reverse=vci_reverse,
                num_processes=None,
                also_patch=True,
                debug_level=verbose
            )

        # normal extraction
        else:
            if identifier:
                g2g_fasta.extract_id(
                    fasta_file=fasta,
                    identifier=identifier,
                    output_file_name=None,
                    reverse=reverse,
                    complement=complement,
                    raw=raw,
                    debug_level=verbose
                )

                return

            elif database:
                if flag_count >= 1:
                    g2g.exit(
                        "Options --complement, --reverse, --reversecomplement "
                        "cannot be used with Database."
                    )

                if transcripts:
                    g2g_fasta.fasta_extract_transcripts(
                        fasta_file=fasta,
                        database_file_name=database,
                        output=None,
                        raw=raw,
                        debug_level=verbose
                    )

                    return

                elif genes:
                    logger.info("Extracting genes from database...")
                    genes = gtf_db.get_genes_simple(
                        database, debug_level=verbose
                    )

                    logger.debug("Genes extracted, transforming to locations")
                    for gene in genes:
                        r = g2g.Region(
                            gene.seqid,
                            gene.start - 1,
                            gene.end,
                            gene.strand,
                            name=gene.ensembl_id,
                            original_base=1,
                        )
                        all_regions.append(r)

                elif exons:
                    exon_ids = {}
                    transcripts = gtf_db.get_transcripts_simple(database)
                    for i, transcript in enumerate(transcripts):
                        for ensembl_id, exon in transcript.exons.items():
                            if ensembl_id not in exon_ids:
                                r = g2g.Region(
                                    exon.seqid,
                                    exon.start - 1,
                                    exon.end,
                                    exon.strand,
                                    name=exon.ensembl_id,
                                    original_base=1,
                                )
                                all_regions.append(r)
                                exon_ids[ensembl_id] = 1

            g2g_fasta.extract(
                fasta_file=fasta,
                locations=all_regions,
                output_file_name=None,
                reverse=reverse,
                complement=complement,
                raw=raw,
                debug_level=verbose,
            )
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GBedError as e:
        g2g.exit(e.msg)
    except exceptions.G2GRegionError as e:
        g2g.exit(e.msg)
    except exceptions.G2GFastaError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# patch
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Fasta file to extract from",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True)
)
@click.option(
    "-c",
    "--vci",
    help="VCI file to use",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True)
)
@click.option(
    "-o",
    "--output",
    help="Name of output file",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True)
)
@click.option(
    "-p",
    "--num-processes",
    help="The number of processes to use, defaults to the number of cores",
    type=int
)
@click.option(
    "-b",
    "--bed",
    help="BED file name",
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "--region",
    help="Region to extract in chromosome:start-end format",
    type=str
)
@click.option(
    "-r",
    "--reverse",
    default=False,
    help="Reverse the direction of the VCI file",
    is_flag=True
)
@click.option(
    "--bgzip",
    default=False,
    help="Compress and index output",
    is_flag=True
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def patch(ctx, input, vci, output, num_processes, bed, region, reverse,
          bgzip, verbose):
    """
    Patch SNPs onto the reference sequence.
    """
    if not output and bgzip:
        g2g.exit("--bgzip can only be used with the -o option")

    all_locations = None

    if region and bed:
        g2g.exit("Please use either a location or a BED file.")
    elif region and not bed:
        all_locations = [g2g.parse_region(region)]
    elif not region and bed:
        bed_file = g2g_bed.BED(bed)
        all_locations = []
        for bed_rec in bed_file:
            if bed_file.current_line_is_bed:
                strand = bed_rec.strand if bed_rec.strand else "+"
                all_locations.append(
                    g2g.Region(
                        bed_rec.chrom,
                        bed_rec.start,
                        bed_rec.end,
                        strand,
                        name=bed_rec.name,
                    )
                )

    try:
        g2g_fasta_patch.process(
            filename_fasta=input,
            filename_vci=vci,
            regions=all_locations,
            filename_output=output,
            bgzip=bgzip,
            reverse=reverse,
            num_processes=num_processes,
            debug_level=verbose
        )
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GVCFError as e:
        g2g.exit(e.msg)
    except exceptions.G2GError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# transform
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Fasta file to extract from",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-c",
    "--vci",
    help="VCI file to use",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-o",
    "--output",
    help="Name of output file",
    type=click.Path(
        exists=False, dir_okay=False, writable=True, resolve_path=True
    )
)
@click.option(
    "-p",
    "--num-processes",
    help="The number of processes to use, defaults to the number of cores",
    type=int
)
@click.option(
    "-b",
    "--bed",
    help="BED file name",
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "--region",
    help="Region to extract in chromosome:start-end format",
    type=str
)
@click.option(
    "-r",
    "--reverse",
    default=False,
    help="Reverse the direction of the VCI file",
    is_flag=True
)
@click.option(
    "--bgzip",
    default=False,
    help="Compress and index output",
    is_flag=True
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def transform(ctx, input, vci, output, num_processes, bed, region, reverse,
              bgzip, verbose):
    """
    Incorporate indels onto the input sequence.
    """
    try:
        if not output and bgzip:
            g2g.exit("--bgzip can only be used with the -o option")

        all_locations = None

        if region and bed:
            g2g.exit("Please use either a location or a BED file.")
        elif region and not bed:
            all_locations = [g2g.parse_region(region)]
        elif not region and bed:
            bed_file = g2g_bed.BED(bed)
            all_locations = []
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else "+"
                    all_locations.append(
                        g2g.Region(
                            bed_rec.chrom, bed_rec.start, bed_rec.end, strand
                        )
                    )

        g2g_fasta_transform.process(
            fasta_file_name_in=input,
            vci_file_name=vci,
            regions=all_locations,
            fasta_file_name_out=output,
            bgzip=bgzip,
            reverse=reverse,
            num_processes=num_processes,
            also_patch=False,
            debug_level=verbose,
        )
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GRegionError as e:
        g2g.exit(e.msg)
    except exceptions.G2GFastaError as e:
        g2g.exit(e.msg)
    except exceptions.G2GError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# parse_region
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-r",
    "--region",
    help="Region to extract in chromosome:start-end format",
    required=True,
    type=str
)
@click.option(
    "-b",
    "--base",
    default="1",
    help="Base coordinate (0 or 1)",
    type=click.Choice(["0", "1"], case_sensitive=False)
)
@click.option(
    "--name",
    help="Region to extract in chromosome:start-end format",
    type=str
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def parse_region(ctx, region, base, name, verbose):
    """
    Parse a region from the command line.  Useful for verifying regions.
    """
    try:
        logger = g2g.get_logger(verbose)

        logger.debug(f"Input: {region}")

        if name:
            logger.debug(f"-n '{name}'")

        base_fixed = 1

        if base:
            logger.debug(f"-b {base}")
            base_fixed = int(base)

        region = g2g.parse_region(region, base=base_fixed, name=name)
        logger.warn(region)
        logger.warn(f"Seq ID: {region.seq_id}")
        logger.warn(f"Start: {region.start}")
        logger.warn(f"_Start: {region._start}")
        logger.warn(f"End: {region.end}")
        logger.warn(f"Name: {region.name}")
        logger.warn(f"Base: {region.original_base}")
        logger.warn(f"Display Start: {region.get_start()}")

    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GBedError as e:
        g2g.exit(e.msg)
    except exceptions.G2GRegionError as e:
        g2g.exit(e.msg)
    except exceptions.G2GFastaError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# fasta_format
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-f",
    "--fasta",
    help="Input Fasta file",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-l",
    "--length",
    default=60,
    help="The length of the line, defaults to 60",
    type=int
)
@click.option(
    "-o",
    "--output",
    help="Output Fasta file",
    required=True,
    type=click.Path(
        exists=False, dir_okay=False, writable=True, resolve_path=True
    )
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def fasta_format(ctx, fasta, length, output, verbose):
    """
    Reformat a Fasta file.
    """
    try:
        g2g_fasta.reformat(
            fasta_file_name=fasta,
            output_file_name=output,
            length=length,
            debug_level=verbose
        )
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GVCFError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# gtf2db
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-i",
    "--input-file",
    help="Input GTF file",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-o",
    "--output-file",
    help="Output database file",
    required=True,
    type=click.Path(
        exists=False, dir_okay=False, writable=True, resolve_path=True
    )
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def gtf2db(ctx, input_file, output_file, verbose):
    """
    Convert a GTF file to a G2G DB file
    """
    try:
        gtf_db.gtf2db(input_file, output_file, verbose)
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GError as e:
        g2g.exit(e.msg)


# #############################################################################
#
# vciquery
#
# #############################################################################
@cli.command()
@click.pass_context
@click.option(
    "-c",
    "--vci",
    help="Input VCI file",
    required=True,
    type=click.Path(
        exists=True, dir_okay=False, readable=True, resolve_path=True
    )
)
@click.option(
    "-r",
    "--region",
    help="Region to extract in chromosome:start-end format",
    required=True,
    type=str
)
@click.option("-v", "--verbose", count=True, help="verbose output")
def vciquery(ctx, vci, region, verbose):
    """
    Query a VCI file.
    """
    try:
        region = g2g.parse_region(region, 1)
        g2g_vci.vci_query(vci, region, debug_level=verbose)
    except KeyboardInterrupt as ki:
        g2g.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g.exit(e.msg)
    except exceptions.G2GError as e:
        g2g.exit(e.msg)

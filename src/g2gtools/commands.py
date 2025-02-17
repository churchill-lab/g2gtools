# standard library imports
import logging
from enum import Enum
from pathlib import Path
from typing import Annotated
import importlib.metadata

# 3rd party library imports
import typer

# local library imports
import g2gtools.bcsam as bcsam
import g2gtools.bed as g2g_bed
import g2gtools.exceptions as exceptions
import g2gtools.fasta as g2g_fasta
import g2gtools.fasta_patch as g2g_fasta_patch
import g2gtools.fasta_transform as g2g_fasta_transform
import g2gtools.g2g_utils as g2g_utils
import g2gtools.gtf as gtf
import g2gtools.gtf_db as gtf_db
import g2gtools.region as g2g_region
import g2gtools.vci as g2g_vci
import g2gtools.vcf2vci as g2g_vcf2vci


class FileFormatEnum(str, Enum):
    BAM = 'BAM'
    SAM = 'SAM'
    GFF = 'GFF'
    GTF = 'GTF'
    BED = 'BED'

app = typer.Typer(help='g2gtools')


def version_callback(value: bool):
    if value:
        from g2gtools import __logo_text__ as logo
        version = importlib.metadata.version('g2gtools')
        typer.echo(f'{logo} {version}')
        raise typer.Exit()


@app.callback()
def common(
    ctx: typer.Context,
    version: bool = typer.Option(None, "--version", callback=version_callback),
):
    pass

# #############################################################################
#
# convert
#
# #############################################################################
@app.command(help='Convert coordinates of BAM|SAM|GTF|GFF|BED files')
def convert(
    input_file: Annotated[Path, typer.Option('--in', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Input file to convert to new coordinates')],
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI file')],
    output_file: Annotated[Path, typer.Option('--out', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')],
    file_format: Annotated[FileFormatEnum, typer.Option('-f', '--file-format',show_default=False,  help='Input file format', case_sensitive=False)],
    reverse: Annotated[bool, typer.Option('--reverse', help='Reverse the direction of the conversion')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0
) -> None:
    """
    Convert coordinates of BAM|SAM|GTF|GFF|BED files
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('convert')

    try:
        input_file = str(input_file) if input_file else None
        vci_file = str(vci_file) if vci_file else None
        output_file = str(output_file) if output_file else None

        if file_format.value.upper() in ['BAM', 'SAM']:
            bcsam.convert_bcsam_file(
                vci_file=vci_file,
                bcsam_file_name_in=input_file,
                bcsam_file_name_out=output_file,
                reverse=reverse
            )
        elif file_format.value.upper() in ['GTF']:
            gtf.convert_gtf_file(
                vci_file=vci_file,
                gtf_file_name_in=input_file,
                gtf_file_name_out=output_file,
                reverse=reverse
            )
        elif file_format.value.upper() in ['BED']:
            g2g_bed.convert_bed_file(
                vci_file=vci_file,
                bed_file_name_in=input_file,
                bed_file_name_out=output_file,
                reverse=reverse
            )
        elif file_format.value.upper() in ['GFF']:
            gtf.convert_gff_file(
                vci_file=vci_file,
                gff_file_name_in=input_file,
                gff_file_name_out=output_file,
                reverse=reverse
            )
        else:
            logger.error(f'Unknown file format: {file_format}')
            raise typer.Exit()
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GChainFileError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBAMError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# convert
#
# #############################################################################
@app.command(help='Convert a region to the new coordinates')
def convert_region(
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI file')],
    region: Annotated[list[str], typer.Option('-r', '--region', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Region to parse in chromosome:start-end format')],
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0
) -> None:
    """
    Convert coordinates of a region
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('convert-region')

    logger.debug(f'reg={region}')

    try:
        regions = g2g_vci.convert_region(vci_file=str(vci_file), reg=region)
        for r in regions:
            logger.warning(f'{r["original"]} -> {r["new"]}')

    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GRegionError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GChainFileError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBAMError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# vcf2vci
#
# #############################################################################
@app.command(help='Create VCI file from VCF file(s)')
def vcf2vci(
    vcf_files: Annotated[list[Path], typer.Option('--vcf', show_default=False, exists=False, dir_okay=False, resolve_path=True, help='VCF files can seperate files by "," or have multiple -vcf')],
    fasta_file: Annotated[Path, typer.Option('--fasta', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Fasta file matching VCF information')],
    strain: Annotated[str, typer.Option('--strain', show_default=False, help='Name of strain/sample (column in VCF file)')],
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')],
    num_processes: Annotated[int, typer.Option('--num-processes', hidden=True)] = None,
    diploid: Annotated[bool, typer.Option('--diploid', help='Create diploid VCI file')] = False,
    keep: Annotated[bool, typer.Option('--keep', help='Keep track of VCF lines that could not be converted to VCI file')] = False,
    passed: Annotated[bool, typer.Option('--pass', help='Use only VCF lines that have a PASS for the filter value')] = False,
    quality: Annotated[bool, typer.Option('--quality', help='Filter on quality, FI=PASS')] = False,
    no_bgzip: Annotated[bool, typer.Option('--no-bgzip', help='DO NOT compress and index output')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0
) -> None:
    """
    Create VCI file from VCF file(s)
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('vcf2vci')

    try:
        # file shortcut: the following command line options are all equal
        # -i abc.vcf -i def.vcf
        # -i abc.vcf,def.vcf
        all_vcf_files: list[str] = []
        for x in vcf_files:
            all_vcf_files.extend(str(x).split(','))

        for i, f in enumerate(all_vcf_files):
            all_vcf_files[i] = g2g_utils.check_file(f, 'r')

        fasta_file = str(fasta_file) if fasta_file else None
        vci_file = str(vci_file) if vci_file else None

        g2g_vcf2vci.process(
            vcf_files=all_vcf_files,
            fasta_file=fasta_file,
            output_file=vci_file,
            strain=strain,
            vcf_keep=keep,
            passed=passed,
            quality=quality,
            diploid=diploid,
            num_processes=num_processes,
            bgzip=(not no_bgzip)
        )
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# extract
#
# #############################################################################
@app.command(help='Extract subsequence from a fasta file given a region')
def extract(
    fasta_file: Annotated[Path, typer.Option('--fasta', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Fasta file to extract from')],
    db_file: Annotated[Path, typer.Option('--db', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Database file name, use with --genes, --transcripts, --exons')] = None,
    out_file: Annotated[Path, typer.Option('--out', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')] = None,
    bed_file: Annotated[Path, typer.Option('-b', '--bed', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='BED file name')] = None,
    genes: Annotated[bool, typer.Option('--genes', help='Extract genes from --database')] = False,
    transcripts: Annotated[bool, typer.Option('--transcripts', help='Extract transcripts from --database')] = False,
    exons: Annotated[bool, typer.Option('--exons', help='Extract exons from --database')] = False,
    identifier: Annotated[str, typer.Option('-id', '--identifier', show_default=False, help='Fasta identifier')] = None,
    region: Annotated[str, typer.Option('-r', '--region', show_default=False, help='Region to extract in chromosome:start-end format')] = None,
    vci_file: Annotated[Path, typer.Option('-c', '--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI File to use')] = None,
    vci_reverse: Annotated[bool, typer.Option('-R', '--vci-reverse', help='Reverse the direction of the VCI file')] = False,
    complement: Annotated[bool, typer.Option('--complement', help='Complement the extracted sequence')] = False,
    reverse: Annotated[bool, typer.Option('--reverse', help='Reverse the extracted sequence')] = False,
    reverse_complement: Annotated[bool, typer.Option('--reverse-complement', help='Reverse-complement the extracted sequence')] = False,
    raw: Annotated[bool, typer.Option('--raw', help='Shows just the extracted sequences')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0
) -> None:
    """
    Extract subsequence from a fasta file given a region

    Note:
        Locations specified on command line are 1-based coordinates.
        Locations specified via BED file are 0-based coordinates.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('extract')

    fasta_file = str(fasta_file) if fasta_file else None
    db_file = str(db_file) if db_file else None
    out_file = str(out_file) if out_file else None
    bed_file = str(bed_file) if bed_file else None
    vci_file = str(vci_file) if vci_file else None

    flag_count = 0
    flag_count += 1 if complement else 0
    flag_count += 1 if reverse else 0
    flag_count += 1 if reverse_complement else 0

    if flag_count > 1:
        logger.error(
            'Please specify only one of: --complement, --reverse, '
            '--reverse-complement',
        )
        raise typer.Exit()

    reverse = reverse or reverse_complement
    complement = complement or reverse_complement

    if db_file:
        # allow only either exon, gene, or transcript extraction, bot all
        output_count = 0
        output_count += 1 if genes else 0
        output_count += 1 if transcripts else 0
        output_count += 1 if exons else 0

        if output_count != 1:
            logger.error(
                'Please specify one of: --exons, --genes, or --transcripts'
            )
            raise typer.Exit()

    location_count = 0
    location_count += 1 if region else 0
    location_count += 1 if bed_file else 0
    location_count += 1 if db_file else 0
    location_count += 1 if identifier else 0

    if location_count > 1:
        logger.error('Please specify only one of: -r, -b, -id, or -db')
        raise typer.Exit()

    # parse the regions we are extracting
    all_regions = []

    try:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta_file_obj = g2g_fasta.FastaFile(fasta_file)

        # get the regions if there are some
        if region:
            # simple region
            logger.debug(f'ARGS.REGION={region}')
            region = g2g_region.parse_region(region, base=1)

            logger.debug(f'--> start = {region.start}')
            logger.debug(f'--> _start = {region._start}')
            logger.debug(f'--> get_start() = {region.get_start()}')

            logger.debug(f'REGION PARSED={region}')
            all_regions.append(region)

            # TEMPORARY OVERRIDE
            # all_regions = [region.parse_region(args.region)]

        elif bed_file:
            # bed file
            bed_file = g2g_bed.BED(bed_file)
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else "+"
                    logger.debug(bed_rec)
                    all_regions.append(
                        g2g_region.Region(
                            bed_rec.chrom,
                            bed_rec.start,
                            bed_rec.end,
                            strand,
                            bed_rec.name,
                            0,
                        )
                    )

        # on the fly
        if vci_file:
            # let's eliminate what isn't supported

            if fasta_file_obj.is_diploid():
                logger.error('Diploid Fasta files are not supported with -v, --vci')
                raise typer.Exit()
            elif identifier:
                logger.error(
                    'Options -v, --vci cannot be used with --id, --identifier'
                )
                raise typer.Exit()
            elif db_file:
                logger.error(
                    'Options -v, --vci cannot be used with --db, --database'
                )
                raise typer.Exit()
            elif flag_count >= 1:
                logger.error(
                    'Options --complement, --reverse, --reversecomplement '
                    'cannot be used with VCI.'
                )
                raise typer.Exit()

            if region:
                # overridden here, fasta.extract expects 1 base
                all_regions = [g2g_region.parse_region(region)]

            g2g_fasta_transform.process(
                filename_fasta=fasta_file,
                filename_vci=vci_file,
                regions=all_regions,
                filename_output=out_file,
                bgzip=False,
                reverse=vci_reverse,
                num_processes=None,
                also_patch=True
            )

        # normal extraction
        else:
            if identifier:
                g2g_fasta.extract_id(
                    fasta_file=fasta_file,
                    identifier=identifier,
                    output_file_name=out_file,
                    reverse=reverse,
                    complement=complement,
                    raw=raw
                )

                return

            elif db_file:
                if flag_count >= 1:
                    logger.error(
                        'Options --complement, --reverse, --reversecomplement '
                        'cannot be used with Database.'
                    )
                    raise typer.Exit()

                logger.warning(f'Input DB File: {db_file}')

                if transcripts:
                    g2g_fasta.fasta_extract_transcripts(
                        fasta_file=fasta_file,
                        database_file_name=db_file,
                        output_file_name=out_file,
                        raw=raw
                    )

                    return

                elif genes:
                    logger.info('Extracting genes from database...')
                    genes = gtf_db.get_genes_simple(db_file)

                    logger.debug('Genes extracted, transforming to locations')
                    for gene in genes:
                        r = g2g_region.Region(
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
                    transcripts = gtf_db.get_transcripts_simple(db_file)
                    for i, transcript in enumerate(transcripts):
                        for ensembl_id, exon in transcript.exons.items():
                            if ensembl_id not in exon_ids:
                                r = g2g_region.Region(
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
                fasta_file=fasta_file,
                locations=all_regions,
                output_file_name=out_file,
                reverse=reverse,
                complement=complement,
                raw=raw
            )
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GRegionError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GFastaError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# patch
#
# #############################################################################
@app.command(help='Patch SNPs onto the reference sequence')
def patch(
    fasta_file: Annotated[Path, typer.Option('--fasta', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Fasta file to extract from')],
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI File to use')],
    out_file: Annotated[Path, typer.Option('--out', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')],
    num_processes: Annotated[int, typer.Option('--num-processes', hidden=True)] = None,
    bed_file: Annotated[Path, typer.Option('--bed', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='BED file name')] = None,
    region: Annotated[str, typer.Option('--region', show_default=False, help='Region to extract in chromosome:start-end format')] = None,
    reverse: Annotated[bool, typer.Option('--reverse', help='Reverse the direction of VCI file')] = False,
    bgzip: Annotated[bool, typer.Option('--bgzip', help='compress and index output')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0
) -> None:
    """
    Patch SNPs onto the reference sequence.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('patch')

    fasta_file = str(fasta_file) if fasta_file else None
    vci_file = str(vci_file) if vci_file else None
    out_file = str(out_file) if out_file else None
    bed_file = str(bed_file) if bed_file else None

    all_locations = None

    if region and bed_file:
        logger.error('Please use either a location or a BED file.')
        raise typer.Exit()
    elif region and not bed_file:
        all_locations = [g2g_region.parse_region(region)]
    elif not region and bed_file:
        bed_file = g2g_bed.BED(bed_file)
        all_locations = []
        for bed_rec in bed_file:
            if bed_file.current_line_is_bed:
                strand = bed_rec.strand if bed_rec.strand else "+"
                all_locations.append(
                    g2g_region.Region(
                        bed_rec.chrom,
                        bed_rec.start,
                        bed_rec.end,
                        strand,
                        name=bed_rec.name,
                    )
                )

    try:
        g2g_fasta_patch.process(
            filename_fasta=fasta_file,
            filename_vci=vci_file,
            regions=all_locations,
            filename_output=out_file,
            bgzip=bgzip,
            reverse=reverse,
            num_processes=num_processes
        )
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GRegionError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GFastaError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# transform
#
# #############################################################################
@app.command(help='Incorporate indels onto the input sequence')
def transform(
    fasta_file: Annotated[Path, typer.Option('--fasta', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Fasta file to extract from')],
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI File to use')],
    out_file: Annotated[Path, typer.Option('--out', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')],
    num_processes: Annotated[int, typer.Option('--num-processes', hidden=True)] = None,
    bed_file: Annotated[Path, typer.Option('--bed', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='BED file name')] = None,
    region: Annotated[str, typer.Option('--region', show_default=False, help='Region to extract in chromosome:start-end format')] = None,
    reverse: Annotated[bool, typer.Option('--reverse', help='Reverse the direction of VCI file')] = False,
    bgzip: Annotated[bool, typer.Option('--bgzip', help='compress and index output')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0,
) -> None:
    """
    Incorporate indels onto the input sequence.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('transform')

    fasta_file = str(fasta_file) if fasta_file else None
    vci_file = str(vci_file) if vci_file else None
    out_file = str(out_file) if out_file else None
    bed_file = str(bed_file) if bed_file else None

    try:
        all_locations = None

        if region and bed_file:
            logger.error('Please use either a location or a BED file.')
            raise typer.Exit()
        elif region and not bed_file:
            all_locations = [g2g_region.parse_region(region)]
        elif not region and bed_file:
            bed_file = g2g_bed.BED(bed_file)
            all_locations = []
            for bed_rec in bed_file:
                if bed_file.current_line_is_bed:
                    strand = bed_rec.strand if bed_rec.strand else "+"
                    all_locations.append(
                        g2g_region.Region(
                            bed_rec.chrom, bed_rec.start, bed_rec.end, strand
                        )
                    )

        g2g_fasta_transform.process(
            filename_fasta=fasta_file,
            filename_vci=vci_file,
            regions=all_locations,
            filename_output=out_file,
            bgzip=bgzip,
            reverse=reverse,
            num_processes=num_processes,
            also_patch=False
        )
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GRegionError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GFastaError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# parse_region
#
# #############################################################################
@app.command(help='Parse a region from the command line.  Useful for verifying regions.')
def parse_region(
    region: Annotated[Path, typer.Option('-r', '--region', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='Region to parse in chromosome:start-end format')],
    base: Annotated[int, typer.Option('-b', '--base', show_default=False, help='Base coordinate (0 or 1)')] = 1,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0,
) -> None:
    """
    Parse a region from the command line.  Useful for verifying regions.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('parse_region')

    try:
        logger.debug(f'Input: {region}')

        if base != 0 or base != 1:
            raise exceptions.G2GValueError('Base must be 0 or 1')

        region = g2g_region.parse_region(region, base=base)
        logger.warning(region)
        logger.warning(f'Seq ID: {region.seq_id}')
        logger.warning(f'Start: {region.start}')
        logger.warning(f'_Start: {region._start}')
        logger.warning(f'End: {region.end}')
        logger.warning(f'Name: {region.name}')
        logger.warning(f'Base: {region.original_base}')
        logger.warning(f'Display Start: {region.get_start()}')

    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GBedError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GRegionError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GFastaError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# fasta_format
#
# #############################################################################
@app.command(help='Reformat a Fasta file')
def fasta_format(
    fasta_file: Annotated[Path, typer.Option('-f', '--fasta', exists=True, dir_okay=False, resolve_path=True, help='Fasta file to extract from')],
    length: Annotated[int, typer.Option('-l', '--length', hidden=True)] = None,
    output_file: Annotated[Path, typer.Option('-o', '--out', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='specify multiple times for more verbose output')] = 0,
) -> None:
    """
    Reformat a Fasta file.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('fasta_format')

    try:
        fasta_file = str(fasta_file) if fasta_file else None
        output_file = str(output_file) if output_file else None

        g2g_fasta.reformat(
            fasta_file_name=fasta_file,
            output_file_name=output_file,
            length=length
        )
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# gtf2db
#
# #############################################################################
@app.command(help='Convert a GTF file to a G2G DB file')
def gtf2db(
    gtf_file: Annotated[Path, typer.Option('--gtf', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='GTF file')],
    db_file: Annotated[Path, typer.Option('--db', show_default=False, exists=False, dir_okay=False, writable=True, resolve_path=True, help='Name of output file')],
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0,
) -> None:
    """
    Convert a GTF file to a G2G DB file
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('gtf2db')

    try:
        gtf_file = str(gtf_file) if gtf_file else None
        db_file = str(db_file) if db_file else None

        gtf_db.gtf2db(gtf_file, db_file)
    except KeyboardInterrupt as ki:
        logger.error(str(ki))
        raise typer.Exit()
    except exceptions.KeyboardInterruptError as e:
        logger.error(str(e))
        raise typer.Exit()
    except exceptions.G2GValueError as e:
        logger.error(e.msg)
        raise typer.Exit()
    except exceptions.G2GVCFError as e:
        logger.error(e.msg)
        raise typer.Exit()


# #############################################################################
#
# vci_query
#
# #############################################################################
@app.command(help='Query a VCI file.')
def vci_query(
    vci_file: Annotated[Path, typer.Option('--vci', show_default=False, exists=True, dir_okay=False, resolve_path=True, help='VCI File to use')],
    region: Annotated[str, typer.Option('-r', '--region', show_default=False, help='Region to extract in chromosome:start-end format')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', show_default=False, count=True, help='specify multiple times for more verbose output')] = 0,
) -> None:
    """
    Query a VCI file.
    """
    logger = g2g_utils.configure_logging('g2gtools', verbose)
    logger.debug('vci_query')

    try:
        vci_file = str(vci_file) if vci_file else None

        region = g2g_region.parse_region(region, 1)
        g2g_vci.vci_query(vci_file, region)
    except KeyboardInterrupt as ki:
        g2g_utils.exit(str(ki))
    except exceptions.KeyboardInterruptError as e:
        g2g_utils.exit(str(e))
    except exceptions.G2GValueError as e:
        g2g_utils.exit(e.msg)
    except exceptions.G2GError as e:
        g2g_utils.exit(e.msg)

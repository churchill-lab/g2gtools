"""
Collection of functions related to BAM and SAM files.

This module provides utilities for working with BAM and SAM files through pysam.
It includes functions for reading, writing, and manipulating alignment data.

Note: pysam uses 0-based coordinates while many genomic formats use 1-based
coordinates.
"""

# standard library imports
from collections import namedtuple
from typing import Any
import importlib.metadata
import os
import pathlib
import re
import time

# 3rd party library imports
import pysam

# local library imports
from g2gtools.exceptions import G2GError
from g2gtools.exceptions import G2GBAMError
from g2gtools.exceptions import G2GCigarFormatError
from g2gtools.exceptions import KeyboardInterruptError
import g2gtools.g2g_utils as g2g_utils
import g2gtools.vci as vci

logger = g2g_utils.get_logger('g2gtools')

FLAG_NONE = 0x0         # base value
FLAG_PAIRED = 0x1       # read paired
FLAG_PROPER_PAIR = 0x2  # read mapped in proper pair
FLAG_UNMAP = 0x4        # read unmapped
FLAG_MUNMAP = 0x8       # mate unmapped
FLAG_REVERSE = 0x10     # read reverse strand
FLAG_MREVERSE = 0x20    # mate reverse strand
FLAG_READ1 = 0x40       # first in pair
FLAG_READ2 = 0x80       # second in pair
FLAG_SECONDARY = 0x100  # not primary alignment
FLAG_QCFAIL = 0x200     # read fails platform/vendor quality checks
FLAG_DUP = 0x400        # read is PCR or optical duplicate
FLAG_SUPP = 0x800       # supplementary alignment

REGEX_CIGAR = re.compile(r'(\d+)([\w=])')
REGEX_CIGAR_LENGTH = re.compile(r'\D')

CIGAR_M = 'M' # Match
CIGAR_I = 'I' # Insertion
CIGAR_D = 'D' # Deletion
CIGAR_N = 'N' # Reference Skip
CIGAR_S = 'S' # Soft Clipping
CIGAR_H = 'H' # Hard Clipping
CIGAR_P = 'P' # Padding
CIGAR_E = '=' # Match (Exact)
CIGAR_X = 'X' # Mismatch

CIGAR_m = 0
CIGAR_i = 1
CIGAR_d = 2
CIGAR_n = 3
CIGAR_s = 4
CIGAR_h = 5
CIGAR_p = 6
CIGAR_e = 7
CIGAR_x = 8

CIGAR_N2C: dict[str | int : str] = {
    0: 'M',  # alignment match (can be a sequence match or mismatch)
    1: 'I',  # insertion to the reference
    2: 'D',  # deletion from the reference
    3: 'N',  # skipped region from the reference
    4: 'S',  # soft clipping (clipped sequences present in SEQ)
    5: 'H',  # hard clipping (clipped sequences NOT present in SEQ)
    6: 'P',  # padding (silent deletion from padded reference)
    7: '=',  # sequence match
    8: 'X',  # sequence mismatch
    '0': 'M',
    '1': 'I',
    '2': 'D',
    '3': 'N',
    '4': 'S',
    '5': 'H',
    '6': 'P',
    '7': '=',
    '8': 'X',
}

CIGAR_C2N: dict[str:int] = {
    'M': 0,
    'I': 1,
    'D': 2,
    'N': 3,
    'S': 4,
    'H': 5,
    'P': 6,
    '=': 7,
    'X': 8,
}

Cigar = namedtuple('Cigar', ['code', 'length', 'start', 'end'])


def open_alignment_file(file_name: str) -> pysam.AlignmentFile | None:
    """
    Open a sequence alignment file with the appropriate mode based on file
    extension.

    This function detects the file type (BAM, SAM, or CRAM) based on the file
    extension and opens it with the correct pysam mode. It provides a convenient
    way to handle different alignment file formats with a single interface.

    Args:
        file_name: Path to the alignment file to open. Supported formats are:
                  - BAM (.bam): Binary Alignment Map format
                  - SAM (.sam): Sequence Alignment Map text format
                  - CRAM (.cram): Compressed alignment format

    Returns:
        A pysam.AlignmentFile object if the file extension is recognized,
        or None if the file extension is not supported.

    Note:
        The returned file should be closed by the caller when no longer needed.
    """
    # Check if the file exists
    if not os.path.exists(file_name):
        return None

    # Define mapping of file extensions to pysam modes
    mode_map: dict[str, str] = {
        '.bam':  'rb',  # read binary
        '.sam':  'r',   # read text
        '.cram': 'rc',  # read CRAM
    }

    # Extract the file extension and convert to lowercase once
    file_extension = pathlib.Path(file_name).suffix.lower()

    # Get the appropriate mode or return None for unsupported extensions
    mode = mode_map.get(file_extension)
    if mode is None:
        return None

    try:
        return pysam.AlignmentFile(file_name, mode)
    except Exception:
        # Handle potential pysam errors
        return None


def convert_bcsam_file(
    vci_file: str | vci.VCIFile,
    bcsam_file_name_in: str,
    bcsam_file_name_out: str,
    reverse: bool = False,
) -> None:
    """
    Convert genome coordinates (in BAM/CRAM/SAM format) between assemblies.

    Args
        vci_file: Name of the VCI file or a VCIFile object.
        bcsam_file_name_in: Input BAM/CRAM/SAM file to convert.
        bcsam_file_name_out: Name of output BAM/CRAM/SAM file.
        reverse: True to process VCI in reverse.
    """
    start = time.time()
    try:
        if not isinstance(vci_file, vci.VCIFile):
            vci_file = g2g_utils.check_file(vci_file)
            logger.info('Parsing VCI File')
            vci_file = vci.VCIFile(vci_file, reverse=reverse)
            vci_file.parse(reverse=reverse)
            logger.info('VCI file parsed')

        bcsam_file_name_in = g2g_utils.check_file(bcsam_file_name_in)
        bcsam_file_name_out = g2g_utils.check_file(bcsam_file_name_out, 'w')
        unmapped_file_name = f'{bcsam_file_name_out}.unmapped'
        alignment_file_in = open_alignment_file(bcsam_file_name_in)

        logger.warning(f'Input VCI File: {vci_file.file_name}')
        logger.warning(f'Input BAM/CRAM/SAM File: {bcsam_file_name_in}')
        logger.warning(f'Output BAM/CRAM/SAM File: {bcsam_file_name_out}')
        logger.warning(f'Output UNMAPPED File: {unmapped_file_name}')
        logger.info('Converting BAM file')

        new_header = alignment_file_in.header.to_dict()

        # replace 'HD'
        new_header['HD'] = {'VN': 1.0, 'SO': 'coordinate'}

        # replace SQ
        # Filter the records in alignment_file_in.header['SQ'] to keep only
        # those that are in the vci file
        filtered_header_sq = [
            record
            for record in alignment_file_in.header['SQ']
            if record['SN'] in vci_file.contigs
        ]

        new_header['SQ'] = filtered_header_sq

        if 'PG' not in new_header:
            new_header['PG'] = []

        new_header['PG'].append(
            {'ID': 'g2gtools', 'VN': importlib.metadata.version('g2gtools')}
        )

        if 'CO' not in new_header:
            new_header['CO'] = []

        new_header['CO'].append(f'Original file: {bcsam_file_name_in}')
        new_header['CO'].append(f'VCI File: {vci_file.filename}')

        temp_dir, temp_file_name = os.path.split(bcsam_file_name_out)
        parts = temp_file_name.split('.')
        ext = parts[-1]

        if ext.lower() == 'bam':
            new_file = pysam.AlignmentFile(
                bcsam_file_name_out, 'wb', header=new_header
            )
            new_file_unmapped = pysam.AlignmentFile(
                unmapped_file_name, 'wb', template=alignment_file_in
            )
        elif ext.lower() == 'sam':
            new_file = pysam.AlignmentFile(
                bcsam_file_name_out, 'wh', header=new_header
            )
            new_file_unmapped = pysam.AlignmentFile(
                unmapped_file_name, 'wh', template=alignment_file_in
            )
        else:
            raise G2GBAMError(
                'Unable to create file based upon file extension'
            )

        # create map from contig name to new ID index.
        # If build using input BAM, IDs are out of sync and it contains contigs not in the VCI file.
        name_to_id = {}
        for ref_name in vci_file.contigs:
            name_to_id[ref_name] = new_file.get_tid(ref_name)

        total = 0
        total_unmapped = 0
        total_fail_qc = 0

        map_statistics = {
            'total': 0,
            'fail_cannot_map': 0,
            'success_simple': 0,
            'success_complex': 0,
        }

        map_statistics_pair = {
            'total': 0,
            'fail_cannot_map': 0,
            'success_1_fail_2_simple': 0,
            'success_1_fail_2_complex': 0,
            'success_1_simple_2_fail': 0,
            'success_1_simple_2_simple': 0,
            'success_1_simple_2_complex': 0,
            'success_1_complex_2_fail': 0,
            'success_1_complex_2_simple': 0,
            'success_1_complex_2_complex': 0,
        }

        try:
            while True:
                if total and total % 10000 == 0:
                    status_success = 0
                    status_failed = 0

                    for k, v in map_statistics_pair.items():
                        if k.startswith('success'):
                            status_success += v
                        elif k.startswith('fail'):
                            status_failed += v

                    logger.info(
                        f'Processed {total:,} reads, '
                        f'{status_success:,} successful, '
                        f'{status_failed:,} failed'
                    )

                alignment = next(alignment_file_in)
                alignment_new = pysam.AlignedSegment()
                read_chr = alignment_file_in.get_reference_name(
                    alignment.reference_id
                )
                total += 1

                logger.debug('~' * 80)
                logger.debug(
                    f'Converting {alignment.query_name} {read_chr} '
                    f'{alignment.reference_start} {alignment.cigarstring}'
                )

                if alignment.is_qcfail:
                    logger.debug('\tFail due to qc of old alignment')
                    new_file_unmapped.write(alignment)
                    total_fail_qc += 1
                    continue

                if alignment.is_unmapped:
                    logger.debug('\tFail due to unmapped old alignment')
                    new_file_unmapped.write(alignment)
                    total_unmapped += 1
                    continue

                if not alignment.is_paired:
                    logger.debug('SINGLE END ALIGNMENT')
                    map_statistics['total'] += 1

                    alignment_new.query_sequence = alignment.query_sequence
                    alignment_new.flag = FLAG_NONE
                    alignment_new.mapping_quality = alignment.mapping_quality
                    alignment_new.query_name = alignment.query_name
                    alignment_new.query_qualities = alignment.query_qualities
                    alignment_new.set_tags(alignment.get_tags())

                    read_start = alignment.reference_start
                    read_end = alignment.reference_end

                    mappings = vci_file.find_mappings(
                        read_chr, read_start, read_end
                    )

                    # unmapped
                    if mappings is None:
                        logger.debug('\tFail due to no mappings')
                        new_file_unmapped.write(alignment)
                        map_statistics['fail_cannot_map'] += 1

                    elif len(mappings) == 1:
                        if alignment.is_reverse:
                            alignment_new.flag |= FLAG_REVERSE

                        alignment_new.reference_id = name_to_id[
                            mappings[0].to_chr
                        ]
                        alignment_new.reference_start = mappings[0].to_start
                        alignment_new.cigartuples = alignment.cigartuples
                        new_file.write(alignment_new)

                        logger.debug(
                            f'\tSuccess (simple): {alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        map_statistics['success_simple'] += 1

                    else:
                        logger.debug(f'MAPPINGS: {len(mappings)}')
                        for m in mappings:
                            logger.debug(f'> {m}')

                        if alignment.is_reverse:
                            alignment_new.flag |= FLAG_REVERSE

                        alignment_new.reference_id = name_to_id[
                            mappings[0].to_chr
                        ]
                        alignment_new.reference_start = mappings[0].to_start
                        alignment_new.cigartuples = convert_cigar(
                            alignment.cigartuples,
                            read_chr,
                            vci_file,
                            alignment.query_sequence,
                            alignment.reference_start,
                        )
                        new_file.write(alignment_new)

                        logger.debug(
                            f'\tSuccess (complex): {alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        map_statistics['success_complex'] += 1

                else:
                    logger.debug('PAIRED END ALIGNMENT')
                    map_statistics_pair['total'] += 1

                    alignment_new.query_sequence = alignment.query_sequence
                    alignment_new.flag = FLAG_PAIRED
                    alignment_new.mapping_quality = alignment.mapping_quality
                    alignment_new.query_name = alignment.query_name
                    alignment_new.query_qualities = alignment.query_qualities
                    alignment_new.set_tags(alignment.get_tags())

                    if alignment.is_read1:
                        alignment_new.flag |= FLAG_READ1
                    if alignment.is_read2:
                        alignment_new.flag |= FLAG_READ2
                    if alignment.is_reverse:
                        alignment_new.flag |= FLAG_REVERSE
                    if alignment.mate_is_reverse:
                        alignment_new.flag |= FLAG_MREVERSE
                    if alignment.is_proper_pair:
                        alignment_new.flag |= FLAG_PROPER_PAIR

                    read1_chr = alignment_file_in.get_reference_name(
                        alignment.reference_id
                    )
                    read1_start = alignment.reference_start
                    read1_end = alignment.reference_end
                    read1_mappings = vci_file.find_mappings(
                        read1_chr, read1_start, read1_end
                    )

                    read2_mappings = None

                    if alignment.mate_is_unmapped:
                        alignment_new.flag |= FLAG_MUNMAP
                    else:
                        read2_chr = alignment_file_in.get_reference_name(
                            alignment.next_reference_id
                        )
                        read2_start = alignment.next_reference_start
                        read2_end = read2_start + 1

                        try:
                            read2_mappings = vci_file.find_mappings(
                                read2_chr, read2_start, read2_end
                            )
                        except Exception as e:
                            logger.error(e)
                            read2_mappings = None

                    if read1_mappings is None and read2_mappings is None:
                        alignment_new.flag |= FLAG_UNMAP
                        alignment_new.flag |= FLAG_MUNMAP

                        logger.debug('\tFail due to no mappings')
                        new_file_unmapped.write(alignment)
                        map_statistics_pair['fail_cannot_map'] += 1

                    elif (
                        read1_mappings is None
                        and read2_mappings
                        and len(read2_mappings) == 1
                    ):
                        alignment_new.flag |= FLAG_UNMAP
                        alignment_new.reference_start = 0
                        alignment_new.cigarstring = '0M'
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0

                        logger.debug(
                            '\tPair Success (1:fail, 2:simple): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_fail_2_simple'] += 1

                    elif (
                        read1_mappings is None
                        and read2_mappings
                        and len(read2_mappings) > 1
                    ):
                        alignment_new.flag |= FLAG_UNMAP
                        alignment_new.reference_start = 0
                        alignment_new.cigarstring = '0M'
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0

                        logger.debug(
                            '\tPair Success (1:fail, 2:complex): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_fail_2_complex'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) == 1
                        and read2_mappings is None
                    ):
                        alignment_new.flag |= FLAG_MUNMAP
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigarstring = alignment.cigarstring
                        alignment_new.next_reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = 0
                        alignment_new.template_length = 0

                        logger.debug(
                            '\tPair Success (1:simple,2:fail): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_simple_2_fail'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) == 1
                        and read2_mappings
                        and len(read2_mappings) == 1
                    ):
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigarstring = alignment.cigarstring
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0  # CHECK

                        logger.debug(
                            '\tPair Success (1:simple,2:simple): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_simple_2_simple'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) == 1
                        and read2_mappings
                        and len(read2_mappings) > 1
                    ):
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigarstring = alignment.cigarstring
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0  # CHECK

                        logger.debug(
                            '\tPair Success (1:simple,2:complex): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_simple_2_complex'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) > 1
                        and read2_mappings is None
                    ):
                        alignment_new.flag |= FLAG_MUNMAP
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigartuples = convert_cigar(
                            alignment.cigartuples,
                            read_chr,
                            vci_file,
                            alignment.query_sequence,
                            alignment.reference_start,
                        )
                        alignment_new.next_reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = 0
                        alignment_new.template_length = 0  # CHECK

                        logger.debug(
                            '\tPair Success (1:complex,2:fail): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_complex_2_fail'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) > 1
                        and read2_mappings
                        and len(read2_mappings) == 1
                    ):
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigartuples = convert_cigar(
                            alignment.cigartuples,
                            read_chr,
                            vci_file,
                            alignment.query_sequence,
                            alignment.reference_start,
                        )
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0  # CHECK

                        logger.debug(
                            '\tPair Success (1:complex,2:simple): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_complex_2_simple'] += 1

                    elif (
                        read1_mappings
                        and len(read1_mappings) > 1
                        and read2_mappings
                        and len(read2_mappings) > 1
                    ):
                        alignment_new.reference_id = name_to_id[
                            read1_mappings[0].to_chr
                        ]
                        alignment_new.reference_start = read1_mappings[
                            0
                        ].to_start
                        alignment_new.cigartuples = convert_cigar(
                            alignment.cigartuples,
                            read_chr,
                            vci_file,
                            alignment.query_sequence,
                            alignment.reference_start,
                        )
                        alignment_new.next_reference_id = name_to_id[
                            read2_mappings[0].to_chr
                        ]
                        alignment_new.next_reference_start = read2_mappings[
                            0
                        ].to_start
                        alignment_new.template_length = 0  # CHECK

                        logger.debug(
                            '\tPair Success (1:complex,2:complex): '
                            f'{alignment_new.reference_start} '
                            f'{alignment_new.cigarstring}'
                        )
                        new_file.write(alignment_new)
                        map_statistics_pair['success_1_complex_2_complex'] += 1

                    else:
                        raise G2GBAMError(
                            'Unknown BAM/CRAM/SAM conversion/parse situation'
                        )

        except StopIteration:
            logger.info('All reads processed')

        logger.warning(f'Processed {total:,} entries')
        logger.warning(f'{total_unmapped:,} unmapped entries')
        logger.warning(f'{total_fail_qc:,} failed QC')

        if map_statistics['total'] > 0:
            n = (
                map_statistics['success_simple']
                + map_statistics['success_complex']
            )
            logger.info('')
            logger.info('Mapping Summary Single End')
            logger.info(f"  {map_statistics['total']:>10} TOTAL ENTRIES")
            logger.info('')
            logger.info(f'  {n:>10} TOTAL SUCCESS')
            logger.info(f"  {map_statistics['success_simple']:>10} Simple")
            logger.info(f"  {map_statistics['success_complex']:>10} Complex")
            logger.info('')
            logger.info(
                f"  {map_statistics['fail_cannot_map']:>10} TOTAL FAILURES"
            )
            logger.info(
                f"  {map_statistics['fail_cannot_map']:>10} Cannot Map"
            )

        if map_statistics_pair['total'] > 0:
            total_success = 0
            for k, v in map_statistics_pair.items():
                if k.startswith('success'):
                    total_success += v

            logger.info('')
            logger.info('Mapping Summary Paired End')
            logger.info(f"  {map_statistics_pair['total']:>10} TOTAL ENTRIES")
            logger.info('')
            logger.info(f'  {total_success:>10} TOTAL SUCCESS')
            logger.info(
                f"  {map_statistics_pair['success_1_fail_2_simple']:>10} "
                'Read 1 Failed, Read 2 Simple'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_fail_2_complex']:>10} "
                'Read 1 Failed, Read 2 Complex'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_simple_2_fail']:>10} "
                'Read 1 Simple, Read 2 Failed'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_simple_2_simple']:>10} "
                'Read 1 Simple, Read 2 Simple'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_simple_2_complex']:>10} "
                'Read 1 Simple, Read 2 Complex'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_complex_2_fail']:>10} "
                'Read 1 Complex, Read 2 Failed'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_complex_2_simple']:>10} "
                'Read 1 Complex, Read 2 Simple'
            )
            logger.info(
                f"  {map_statistics_pair['success_1_complex_2_complex']:>10} "
                'Read 1 Complex, Read 2 Complex'
            )
            logger.info('')
            logger.info(
                f"  {map_statistics_pair['fail_cannot_map']:>10} TOTAL FAILURES"
            )
            logger.info(
                f"  {map_statistics_pair['fail_cannot_map']:>10} Cannot Map"
            )
            logger.info('')

        logger.info('BAM/SAM/CRAM File converted')
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')


def cigarlist_to_cigarstring(cigar_list: Cigar | list[tuple[int, int]]) -> str:
    """
    Convert a CIGAR list representation to a CIGAR string.

    Transforms a list of CIGAR operations (either as Cigar objects or tuples)
    into the standard string representation used in SAM/BAM files.

    Examples:
        >>> cigarlist_to_cigarstring([(0, 10), (1, 1), (0, 75), (2, 2), (0, 20)])
        '10M1I75M2D20M'

    Args:
        cigar_list: A CIGAR representation, either:
                   - A Cigar object containing CigarOp elements
                   - A list of tuples where each tuple is (operation_code, length)
                     Operation codes follow the BAM format:
                     0=M (match/mismatch), 1=I (insertion), 2=D (deletion), etc.

    Returns:
        A string representation of the CIGAR in standard SAM/BAM format.

    Raises:
        G2GCigarFormatError: If any CIGAR operation has an invalid code.
        TypeError: If the input is neither a Cigar object nor a list of tuples.
    """
    # handle empty input
    if not cigar_list:
        return ''

    # process based on input type
    if isinstance(cigar_list, Cigar):
        # process Cigar object
        try:
            return ''.join(f'{op.length}{op.code}' for op in cigar_list)
        except (AttributeError, KeyError) as e:
            # identify the problematic operation if possible
            for i, op in enumerate(cigar_list):
                if not hasattr(op, 'code') or not hasattr(op, 'length'):
                    raise G2GCigarFormatError(
                        f'Invalid Cigar operation at position {i}: {op}'
                    )
            # if we can't identify the specific issue, raise a general error
            raise G2GCigarFormatError(f'Invalid Cigar object: {e}')

    elif isinstance(cigar_list, list):
        # process list of tuples
        try:
            return ''.join(
                f'{length}{CIGAR_N2C[code]}' for code, length in cigar_list
            )
        except (IndexError, KeyError, ValueError) as e:
            # find the problematic tuple
            for i, item in enumerate(cigar_list):
                if not isinstance(item, tuple) or len(item) != 2:
                    raise G2GCigarFormatError(
                        f'Invalid tuple at position {i}: {item}, '
                        'expected (code, length)'
                    )
                code, length = item
                if code not in CIGAR_N2C:
                    raise G2GCigarFormatError(
                        f'Invalid CIGAR operation code at position {i}: {code}'
                    )
            # if we can't identify the specific issue, raise a general error
            raise G2GCigarFormatError(f'Error processing CIGAR list: {e}')

    else:
        # handle invalid input type
        raise TypeError(
            f'Expected Cigar object or list of tuples, got {type(cigar_list).__name__}'
        )


def cigar_to_string(cigar: list[tuple[int, int]]) -> str:
    """
    Convert a list of tuples into a cigar string.

    Example:
             [ (0, 10), (1, 1), (0, 75), (2, 2), (0, 20) ]
        =>       10M      1I      75M      2D      20M
        =>                10M1I75M2D20M

    Args:
        cigar: The cigar string as a list of tuples.

    Returns:
        The cigar string.
    """
    return ''.join(f"{length}{CIGAR_N2C.get(op, '?')}" for op, length in cigar)



def cigar_convert(
    cigar: list[Any], chromosome: str, vci_file: vci.VCIFile, position: int = 0
) -> list[Cigar]:
    """
    Convert each CIGAR element to new mappings and construct an array on NEW
    cigar elements.

    Example:
        Depending on the Intervals in the VCI file, let's say we have the
        following CIGAR string: 35M49N65M

        This could get converted into:

        35M ==> 4M150D31M
        49N ==> -1N (remember, surrounding M's are used to find the length
                    of N which is done on next pass)
        65M ==> 65M

        First pass yields: 35M49N65M => 4M150D31M-1N65M

    Args:
        cigar: The cigar as a list with alternating strings and integers.
        chromosome: The chromosome.
        vci_file: The VCI file as a vci.VCIFile object.
        position: The cigar position.

    Returns:
        A list of Cigar objects representing the cigar.
    """
    cigar_new: list[Cigar] = []
    current_pos = position
    cigar_no = 0

    for c in cigar:
        cigar_no += 1
        increment = c[1]

        if c[0] == CIGAR_m:
            new_mappings = vci_file.find_mappings(
                chromosome, current_pos, current_pos + c[1]
            )

            if not new_mappings:
                # logger.debug('Mappings: None')
                cigar_new.append(Cigar(CIGAR_S, c[1], 0, 0))
            elif len(new_mappings) == 1:
                # logger.debug(f'Mappings: Easy: {new_mappings[0]}'))
                cigar_new.append(
                    Cigar(
                        CIGAR_M,
                        new_mappings[0].to_end - new_mappings[0].to_start,
                        new_mappings[0].to_start,
                        new_mappings[0].to_end,
                    )
                )
            else:
                # multiple maps, not so easy
                last = None

                for m in new_mappings:
                    # logger.debug(f'Mappings: Multiple: {m}')
                    if not last:
                        last = m

                        if current_pos < m.from_start:
                            # special case of first match not in interval,
                            # handle accordingly
                            # logger.debug(
                            #     f'Adding 'S', because {current_pos} < '
                            #     f'{m.from_start}'
                            # )
                            cigar_new.append(
                                Cigar(
                                    CIGAR_S, m.from_start - current_pos, 0, 0
                                )
                            )
                    else:
                        if m.from_start != last.from_end:
                            # logger.debug(
                            #     'Adding 'M' and 'I', because '
                            #     f'{m.from_start} != {last.from_end}'
                            # )
                            cigar_new.append(
                                Cigar(
                                    CIGAR_M,
                                    last.to_end - last.to_start,
                                    last.to_start,
                                    last.to_end,
                                )
                            )
                            cigar_new.append(
                                Cigar(
                                    CIGAR_I,
                                    m.from_start - last.from_end,
                                    last.to_start,
                                    last.to_end,
                                )
                            )
                        elif m.to_start != last.to_end:
                            # logger.debug(
                            #     f"Adding 'M' and 'D', because {m.to_start} "
                            #     f"!= {last.to_end}"
                            # )
                            cigar_new.append(
                                Cigar(
                                    CIGAR_M,
                                    last.to_end - last.to_start,
                                    last.to_start,
                                    last.to_end,
                                )
                            )
                            cigar_new.append(
                                Cigar(CIGAR_D, m.to_start - last.to_end, 0, 0)
                            )

                        last = m

                # logger.debug("Adding 'M'")
                cigar_new.append(
                    Cigar(
                        CIGAR_M,
                        last.to_end - last.to_start,
                        last.to_start,
                        last.to_end,
                    )
                )

        elif c[0] == CIGAR_i:
            # logger.debug("Adding 'I' and 'D'")
            cigar_new.append(Cigar(CIGAR_I, c[1], 0, 0))
            cigar_new.append(Cigar(CIGAR_D, -1, 0, 0))
            increment = 0
        elif c[0] == CIGAR_d:
            # logger.debug("Adding 'D'")
            cigar_new.append(Cigar(CIGAR_D, -1, 0, 0))
        elif c[0] == CIGAR_n:
            # logger.debug("Adding 'N'")
            cigar_new.append(Cigar(CIGAR_N, -1, 0, 0))
        elif c[0] in [CIGAR_s, CIGAR_h]:
            # logger.debug(f"Adding '{CIGAR_N2C[c[0]]}'")
            cigar_new.append(Cigar(CIGAR_N2C[c[0]], c[1], 0, 0))
        else:
            # other
            # logger.debug(
            #     f"OTHER CODE '{CIGAR_N2C[c[0]]}' found, looking at "
            #     f"{c} at {current_pos}"
            # )
            raise G2GCigarFormatError(
                f'ERROR: Not handling the values in this cigar string: {cigar}'
            )

        current_pos += increment

        # logger.debug(f"Current CIGAR: {cigar_new}")

    return cigar_new


def cigar_combine_consecutive(cigar: list[Cigar]) -> list[Cigar]:
    """
    Combine consecutive CIGAR elements with the same operation code.

    This function iteratively scans through a list of Cigar objects and combines
    adjacent elements that have the same operation code (e.g., M, I, D, S).
    For example, 2M followed by 3M would be combined into 5M, or 10N followed
    by 5N would become 15N.

    The function continues combining elements until no more consecutive elements
    with the same code exist in the list.

    Args:
        cigar: A list of Cigar objects representing a CIGAR string.
               Each Cigar object contains an operation code, length, and
               genomic coordinates.

    Returns:
        A new list of Cigar objects with consecutive identical operations combined.
        The resulting list will have the same or fewer elements than the input list.

    Example:
        >>> # Assuming Cigar objects with code, length, start, end
        >>> cigar = [
        ...     Cigar('M', 10, 100, 110),
        ...     Cigar('M', 5, 110, 115),
        ...     Cigar('I', 2, 115, 115),
        ...     Cigar('M', 20, 115, 135)
        ... ]
        >>> result = cigar_combine_consecutive(cigar)
        >>> # Result would be:
        >>> # [Cigar('M', 15, 100, 115), Cigar('I', 2, 115, 115), Cigar('M', 20, 115, 135)]
    """
    done = False
    while not done:
        done = True
        i = 0
        # Find the first pair of consecutive elements with the same code
        for i in range(0, len(cigar) - 1):
            # logger.debug(f"{i} = {cigar[i]}")
            # logger.debug(f"{i + 1} = {cigar[i + 1]}")
            if cigar[i].code == cigar[i + 1].code:
                done = False
                break

        # If we found a pair to combine
        if not done:
            cigar_temp: list[Cigar] = []

            # Keep all elements before the pair
            cigar_temp.extend(cigar[:i])

            # Get the two consecutive elements
            cm1 = cigar[i]
            cm2 = cigar[i + 1]

            # Create a new combined element
            # The new element has:
            # - Same code as the original elements
            # - Combined length (sum of both lengths)
            # - Start position from the first element
            # - End position from the second element
            cm_new = Cigar(
                cm1.code, cm1.length + cm2.length, cm1.start, cm2.end
            )

            # Add the combined element
            cigar_temp.append(cm_new)

            # logger.debug(
            #     f"Found consecutive elements {cm1} and {cm2}, combined "
            #     f"into {cm_new}"
            # )

            # Keep all elements after the pair
            cigar_temp.extend(cigar[i + 2 :])

            # Update the cigar list for the next iteration
            cigar = cigar_temp

    return cigar


def cigar_fix_pre_and_post_m(cigar: list[Cigar]) -> list[Cigar]:
    """
    Convert elements before the first M and after the last M to soft clips (S).

    This function handles two specific cases in CIGAR string processing:
    1. Pre-M fix: Converts all elements before the first match (M) to a single soft clip (S)
    2. Post-M fix: Converts all elements after the last match (M) to a single soft clip (S)

    This is necessary because elements like insertions (I) or hard clips (H) that appear
    before the first alignment match or after the last match should be represented as
    soft clips in the final CIGAR string.

    Args:
        cigar: A list of Cigar objects representing a CIGAR string.

    Returns:
        A modified list of Cigar objects with proper soft clips at the beginning and end.

    Example:
        >>> # Input: [I(2), H(3), M(10), D(5), M(15), I(3)]
        >>> # Output: [S(5), M(10), D(5), M(15), S(3)]
        >>> # (where the first I(2)+H(3) becomes S(5), and the last I(3) becomes S(3))
    """
    # Pre-M to S fix: Convert all elements before the first M to a single S
    i = 0
    for i in range(0, len(cigar)):
        if cigar[i].code == CIGAR_M:
            break

    if i != 0:  # If there are elements before the first M
        first_m = i
        length = 0
        # Calculate the total length of elements to convert to soft clips
        for i in range(0, first_m):
            if cigar[i].code in [CIGAR_I, CIGAR_S, CIGAR_H]:
                length += cigar[i].length

        # Create a new CIGAR with a soft clip followed by all elements after the first M
        temp_cigar: list[Cigar] = [Cigar(CIGAR_S, length, 0, 0)]
        temp_cigar.extend(cigar[first_m:])
        cigar = temp_cigar

    # Post-M to S fix: Convert all elements after the last M to a single S
    for i in reversed(range(0, len(cigar))):
        if cigar[i].code == CIGAR_M:
            break

    if i > 0 and i != (
        len(cigar) - 1
    ):  # If there are elements after the last M
        last_m = i
        length = 0
        # Calculate the total length of elements to convert to soft clips
        for i in range(last_m + 1, len(cigar)):
            if cigar[i].code in [CIGAR_M, CIGAR_I, CIGAR_S]:
                length += cigar[i].length

        # Create a new CIGAR with all elements up to the last M followed by a soft clip
        temp_cigar: list[Cigar] = []
        temp_cigar.extend(cigar[: last_m + 1])
        temp_cigar.append(Cigar(CIGAR_S, length, 0, 0))
        cigar = temp_cigar

    return cigar


def cigar_remove_softs_between_m(cigar: list[Cigar]) -> list[Cigar]:
    """
    Remove soft clips (S) that are surrounded by matches (M) in a CIGAR string.

    This function iteratively scans the CIGAR string for soft clip elements that
    are directly preceded and followed by match elements. When found, these soft
    clips are removed from the CIGAR string.

    This is typically done because soft clips between matches are not biologically
    meaningful and likely represent artifacts in the alignment process.

    Args:
        cigar: A list of Cigar objects representing a CIGAR string.

    Returns:
        A modified list of Cigar objects with soft clips between matches removed.

    Example:
        >>> # Input: [M(10), S(5), M(20)]
        >>> # Output: [M(10), M(20)]
    """
    done = False
    while not done:
        done = True
        i = 0

        # Find a soft clip that's not at the beginning or end
        for i in range(1, len(cigar) - 1):
            if cigar[i].code == CIGAR_S:
                done = False
                break

        if done:
            break

        # Look for M elements before and after the soft clip
        before = None
        after = None

        for x in reversed(range(i)):
            if cigar[x].code == CIGAR_M:
                before = cigar[x]
                break

        for x in range(i + 1, len(cigar)):
            if cigar[x].code == CIGAR_M:
                after = cigar[x]
                break

        # If the soft clip is surrounded by matches, remove it
        if before and after:
            # logger.debug("Found 'S' between 'M' so removing 'S'")
            cigar_temp: list[Cigar] = []
            cigar_temp.extend(cigar[:i])
            cigar_temp.extend(cigar[i + 1 :])
            cigar = cigar_temp
            # logger.debug(cigar)
        else:
            done = True

    return cigar


def cigar_fix_lengths(cigar: list[Cigar], sequence: str) -> list[Cigar]:
    """
    Fix CIGAR elements with undefined lengths (-1) and remove zero-length elements.

    This function handles two main tasks:
    1. Assigns proper lengths to CIGAR elements marked with -1 length
       (typically deletions or skipped regions between matches)
    2. Removes any CIGAR elements that end up with zero length

    For elements with -1 length, the function calculates the appropriate length
    based on the genomic coordinates of surrounding match (M) elements. If no
    surrounding matches are found, the element is converted to a soft clip.

    Args:
        cigar: A list of Cigar objects representing a CIGAR string.
        sequence: The read sequence that corresponds to the CIGAR string.
               Used to calculate proper lengths for soft clips.

    Returns:
        A modified list of Cigar objects with all lengths properly defined
        and zero-length elements removed.

    Example:
        >>> # Input with -1 length N between two M elements:
        >>> # [M(35), N(-1), M(65)]
        >>> # Output (assuming the genomic distance is 49):
        >>> # [M(35), N(49), M(65)]
    """
    done = False
    while not done:
        done = True

        # Find first element with undefined length (-1)
        i = 0
        for i, cm in enumerate(cigar):
            if cm.length == -1:
                break

        if i == len(cigar):  # No elements with -1 length found
            done = True
            break

        # logger.debug(f"Found '{cm.code}' at {i}: {cm}")

        # Look for M elements before and after the -1 length element
        before = None
        after = None

        for x in reversed(range(i)):
            if cigar[x].code == CIGAR_M:
                before = cigar[x]
                break

        for x in range(i + 1, len(cigar)):
            if cigar[x].code == CIGAR_M:
                after = cigar[x]
                break

        # Special case handling: check if all remaining elements have undefined length
        a = i
        while a < len(cigar) - 1:
            if cigar[a].length != -1:
                break
            a += 1

        # If we can't determine the length from surrounding matches,
        # convert the rest to a soft clip
        if (
            (a == len(cigar) - 1 and cigar[a].start == -1)
            or not after
            or not before
        ):
            # logger.debug('Found a clip')
            temp_cigar_mappings = cigar[:i]

            # Calculate the length for the soft clip
            temp_total = 0
            for t in temp_cigar_mappings:
                if t.code in [CIGAR_M, CIGAR_I, CIGAR_S]:
                    temp_total += t.length

            # Add a soft clip for the remaining sequence
            temp_cigar_mappings.append(
                Cigar(CIGAR_S, len(sequence) - temp_total, -1, -1)
            )

            cigar = temp_cigar_mappings
            done = True
        else:
            # Calculate the length based on genomic coordinates
            c = cigar[i]
            new_c = Cigar(
                c.code, after.start - before.end, before.end, after.start
            )
            # logger.debug(f"Replacing, old = {c}, new = {new_c}")
            cigar[i] = new_c
            done = False

    # Remove zero-length elements
    # logger.debug("Removing 0 length elements, if any")
    new_cigar: list[Cigar] = []
    for cm in cigar:
        if cm.length == 0:  # Fixed from cm[1] to cm.length
            # logger.debug(f"Removing {cm}")
            continue
        new_cigar.append(cm)

    return new_cigar


def convert_cigar(
    cigar: Cigar | list[tuple[int, int]],
    chromosome: str,
    vci_file: vci.VCIFile,
    sequence: str,
    position: int = 0,
) -> list[tuple[int, int]]:
    """
    Generate a new CIGAR string by converting an alignment between genome assemblies.

    This function transforms a CIGAR string from one genome assembly to another using
    the mapping information in a VCI file. The conversion follows several phases:

    1. Map each CIGAR element to new coordinates
    2. Remove soft clips (S) if surrounded by matches (M)
    3. Fix element lengths for insertions, deletions, and skipped regions
    4. Combine consecutive matching elements
    5. Fix elements before the first match and after the last match
    6. Validate and adjust the final CIGAR string

    Args:
        cigar: CIGAR tuples in the form [(op_code, length), ...] where op_code is an
               integer representing the operation (0=M, 1=I, 2=D, etc.) and length is
               the number of bases affected.
        chromosome: The chromosome name where the alignment is located.
        vci_file: A VCIFile object containing the mapping between genome assemblies.
        sequence: The read sequence that corresponds to the CIGAR string.
        position: The starting position of the alignment on the chromosome (0-based).
                 Defaults to 0.

    Returns:
        A list of CIGAR tuples in the form [(op_code, length), ...] representing
        the converted alignment in the target genome assembly.

    Example:
        >>> vci_file = vci.VCIFile("mouse.vci")
        >>> old_cigar = [(0, 10), (1, 1), (0, 75), (2, 2), (0, 20)]
        >>> new_cigar = convert_cigar(old_cigar, "chr1", vci_file, "ACGTACGTAC...")
    """
    logger = g2g_utils.get_logger()
    old_cigar = cigarlist_to_cigarstring(cigar)
    logger.debug(f'CIGAR CONVERSION : {old_cigar}')

    #
    # PHASE 1: Convert each CIGAR element to new mappings and construct an
    #          array on NEW cigar elements
    #
    logger.debug('CIGAR CONVERSION : PHASE 1 : Converting cigar elements')
    new_cigar = cigar_convert(cigar, chromosome, vci_file, position)
    logger.debug(f'AFTER PHASE 1 : {new_cigar}')

    if len(new_cigar) == 1:
        logger.debug('CIGAR CONVERSION : Skipping to end since only 1 element')
    else:
        #
        # PHASE 2: Remove S if surrounded by M
        #
        logger.debug(
            'CIGAR CONVERSION : PHASE 2 : Remove S if surrounded by M'
        )
        new_cigar = cigar_remove_softs_between_m(new_cigar)
        logger.debug(f'AFTER PHASE 2 : {new_cigar} ')

        #
        # PHASE 3: Fix element lengths
        #
        logger.debug('CIGAR CONVERSION : PHASE 3 : Fix element lengths')
        new_cigar = cigar_fix_lengths(new_cigar, sequence)
        logger.debug(f'AFTER PHASE 3 : {new_cigar} ')

        #
        # PHASE 4: Combine consecutive matching elements
        #
        logger.debug('CIGAR CONVERSION : PHASE 4 : Combining elements')
        new_cigar = cigar_combine_consecutive(new_cigar)
        logger.debug(f'AFTER PHASE 4 : {new_cigar} ')

        #
        # PHASE 5: Fix elements before first M and after last M
        #
        logger.debug('CIGAR CONVERSION : PHASE 5 : Fix pre and post Ms')
        new_cigar = cigar_fix_pre_and_post_m(new_cigar)
        logger.debug(f'AFTER PHASE 5 : {new_cigar} ')

    #
    # Final pass through CIGAR string
    #
    # Test cigar string length
    #
    # SEQ: segment SEQuence. This field can be a '*' when the sequence is not
    #      stored. If not a '*', the length of the sequence must equal the sum
    #      of lengths of M/I/S/=/X operations in CIGAR. An '=' denotes the base
    #      is identical to the reference base. No assumptions can be made on
    #      the letter cases.
    #
    logger.debug('CIGAR CONVERSION : PHASE 6 : Testing length and conversion')
    cigar_seq_length = 0

    # simplify the cigar, throw away the other stuff we used
    simple_cigar: list[tuple[int, int]] = []
    for c in new_cigar:
        simple_cigar.append((CIGAR_C2N[c.code], c.length))
        if c.code in [CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_E, CIGAR_X]:
            cigar_seq_length += c.length

    try:
        if cigar_seq_length != len(sequence):
            logger.debug(
                f'CIGAR SEQ LENGTH={cigar_seq_length} != SEQ_LEN={len(sequence)}'
            )
            # Not equal according to chain file format, add the clipping length
            simple_cigar.append((CIGAR_s, len(sequence) - cigar_seq_length))
    except TypeError:
        # avoids: TypeError: object of type 'NoneType' has no len()
        pass

    if old_cigar != cigar_to_string(simple_cigar):
        logger.debug('old cigar != new cigar')
    else:
        logger.debug('old cigar == new cigar')

    logger.debug(
        f'CIGAR CONVERSION : {old_cigar} ==> {cigar_to_string(simple_cigar)}'
    )
    logger.debug(simple_cigar)

    return simple_cigar


# if __name__ == '__main__':
#    log = g2g_utils.get_logger()
#    cigarstring = '5I3D4M9D3S104M7D2I'
#    cigarlist = cigar_string_to_list(cigarstring)
#    log.debug(cigarstring)
#    log.debug(cigarlist)

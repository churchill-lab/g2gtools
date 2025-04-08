"""
A way to combine vcf files into a more succinct format.
"""
# standard library imports
from dataclasses import dataclass
from collections import OrderedDict
from typing import Any, Callable, Iterator, TypeVar
import copy
import logging
import multiprocessing
import os
import time

# 3rd party library imports
import pysam

# local library imports
from g2gtools.exceptions import G2GError
from g2gtools.exceptions import G2GValueError
from g2gtools.exceptions import G2GVCFError
from g2gtools.exceptions import KeyboardInterruptError
import g2gtools.fasta as fasta
import g2gtools.g2g_utils as g2g_utils
import g2gtools.vcf as vcf
import g2gtools.gtf as gtf

logger = g2g_utils.get_logger('g2gtools')

# Define type variables for VCF readers and records
VCFRecord = TypeVar('VCFRecord')
VCFReader = Iterator[VCFRecord]


@dataclass
class VCFFileInformation:
    """
    Represents information about a VCF file.

    This class stores metadata about a VCF file, including its filename,
    associated discard file, sample name, and sample index.

    Attributes:
        file_name (str): The path to the VCF file.
        discard_file (str | None): Path to a discard file associated with this VCF.
        sample_name (str | None): Name of the sample in the VCF file.
        sample_index (int | None): Index of the sample in a multi-sample VCF file.
    """
    file_name: str
    discard_file: str | None
    sample_name: str | None
    sample_index: int | None

    def __init__(
        self,
        file_name: str,
        discard_file: str | None = None,
    ) -> None:
        """
        Initialize a new VCFFileInformation instance.

        Args:
            file_name: The path to the VCF file.
            discard_file: Path to a discard file associated with this VCF.
                          Defaults to None.

        Note:
            sample_name and sample_index can be set after initialization.
        """
        self.file_name = file_name
        self.discard_file = discard_file
        self.sample_name = None
        self.sample_index = None

    def __str__(self) -> str:
        """
        Return a string representation of the VCF file information.

        Returns:
            str: A string containing the file name and sample name, separated
                 by a comma.
        """
        return f'{self.file_name}, {self.sample_name}'



class VCF2VCIInfo:
    """
    Configuration and state information for VCF to VCI conversion.

    This class stores all the necessary information and settings required for
    converting VCF files to VCI format, including input files, output paths,
    conversion parameters, and statistics tracking.

    Attributes:
        chromosome (str | None): Target chromosome for conversion.
        vcf_files (list[VCFFileInformation]): List of VCF files to process.
        fasta_file (str | None): Path to the reference FASTA file.
        log_level (int | None): Logging level, defaults to logging.INFO.
        gtf_file (str | None): Path to the GTF annotation file.
        prev_next_ref_pos_left (int): Previous reference position for left alignment.
        prev_next_ref_pos_right (int): Previous reference position for right alignment.
        diploid (bool): Whether to process as diploid genome.
        passed (bool): Whether to only include variants that passed filters.
        quality (bool): Whether to include quality information.
        vcf_keep (bool): Whether to keep original VCF files after conversion.
        output_file_left (str | None): Path to the output file for left-aligned variants.
        output_file_right (str | None): Path to the output file for right-aligned variants.
        output_stats_file_left (str | None): Path to the statistics file for left-aligned variants.
        output_stats_file_right (str | None): Path to the statistics file for right-aligned variants.
        stats_left (dict[str, int]): Statistics for left-aligned variants.
        stats_right (dict[str, int]): Statistics for right-aligned variants.
    """

    chromosome: str | None
    vcf_files: list[VCFFileInformation]
    fasta_file: str | None
    log_level: int | None
    gtf_file: str | None
    prev_next_ref_pos_left: int
    prev_next_ref_pos_right: int
    diploid: bool
    passed: bool
    quality: bool
    vcf_keep: bool
    output_file_left: str | None
    output_file_right: str | None
    output_stats_file_left: str | None
    output_stats_file_right: str | None
    stats_left: dict[str, int]
    stats_right: dict[str, int]

    def __init__(self) -> None:
        """
        Initialize a new VCF2VCIInfo instance with default values.

        All attributes are initialized to their default values:
        - File paths are set to None
        - Boolean flags are set to False
        - Containers are initialized as empty
        - Position counters are set to 1
        - Log level defaults to INFO
        """
        self.chromosome = None
        # list of VCFFileInformation
        self.vcf_files = []
        self.fasta_file = None
        self.log_level = logging.INFO
        self.gtf_file = None
        # indel specific
        self.prev_next_ref_pos_right = 1
        self.prev_next_ref_pos_left = 1
        self.diploid = False
        self.passed = False
        self.quality = False
        self.vcf_keep = False
        self.output_file_left = None
        self.output_file_right = None
        self.output_stats_file_left = None
        self.output_stats_file_right = None
        self.stats_left = {}
        self.stats_right = {}


def walk_vcfs_together(
    readers: list[VCFReader | None],
    **kwargs: Any,
) -> Iterator[list[VCFRecord | None]]:
    """
    Simultaneously iterate over two or more VCF readers.

    For each genomic position with a variant, return a list of size equal to the
    number of VCF readers. This list contains the VCF record from readers that have
    this variant, and None for readers that don't have it. Inputs must be
    sorted in the same way and use the same reference.

    Args:
        readers: A list of VCF readers (iterators over VCF records).
            Readers can be None, in which case None will always be returned
            for that position in the result.

        **kwargs: Additional keyword arguments.
            vcf_record_sort_key: If provided, should be a function that takes a
                VCF record and returns a tuple that can be used as a key for
                comparing and sorting VCF records across all readers. This tuple
                defines what it means for two variants to be equal (i.e., whether
                it's only their position or also their allele values), and implicitly
                determines the chromosome ordering since the tuple's 1st element is
                typically the chromosome name (or calculated from it).

    Returns:
        An iterator yielding lists of VCF records (or None). Each list has the same
        length as the input readers list. For each position, the list contains the
        VCF record from each reader that has a variant at that position, or None
        for readers that don't have a variant there.
    """
    # define the default sort key function if not provided
    if 'vcf_record_sort_key' in kwargs:
        get_key: Callable[[VCFRecord], tuple[Any, ...]] = kwargs[
            'vcf_record_sort_key'
        ]
    else:

        def get_key(r: VCFRecord) -> tuple[str, int]:
            return (
                r.contig,
                r.pos,
            )

    # initialize with the first record from each reader
    nexts: list[VCFRecord | None] = []
    for reader in readers:
        try:
            if reader:
                nexts.append(next(reader))
            else:
                nexts.append(None)
        except StopIteration:
            nexts.append(None)

    # start with a minimal key that will be replaced
    min_k: tuple[Any, ...] = (None,)

    # Continue until all readers are exhausted
    while any([r is not None for r in nexts]):
        # map reader indices to their current record's sort key
        next_idx_to_k: dict[int, tuple[Any, ...]] = dict(
            (i, get_key(r)) for i, r in enumerate(nexts) if r is not None
        )

        # find keys that match the previous contig (chromosome)
        keys_with_prev_contig: list[tuple[Any, ...]] = [
            k for k in next_idx_to_k.values() if k[0] == min_k[0]
        ]

        # determine the minimum key (position) to process next
        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)
        else:
            min_k = min(next_idx_to_k.values())

        # find all readers that have a record at this minimum position
        min_k_idxs: set[int] = set(
            [i for i, k in next_idx_to_k.items() if k == min_k]
        )

        # yield a list with records for readers at this position, None for others
        yield [
            nexts[i] if i in min_k_idxs else None for i in range(len(nexts))
        ]

        # advance the readers that were just processed
        for i in min_k_idxs:
            try:
                nexts[i] = next(readers[i])
            except StopIteration:
                nexts[i] = None


def update_stats(stats: dict[str, int], reason: str) -> dict[str, int]:
    """
    Update the statistics to keep track of the processing.

    Args:
        stats: A dictionary of the statistics.
        reason: A string representing the key to the dictionary.

    Returns:
        The statistics.
    """
    if stats:
        stats[reason] = stats.get(reason, 0) + 1
    else:
        stats = {reason: 1}

    return stats


def process_piece(params: VCF2VCIInfo) -> dict[str, Any]:
    """
    Process this 'piece' of the VCF file.

    Args:
        params: The parameters dictating what piece to process.

    Returns:
        The results of the processing.
    """
    logger = g2g_utils.configure_logging('g2gtools', params.log_level)

    stats = {}

    try:
        output_file_left = None
        output_file_right = None
        output_stats_file_left = None
        output_stats_file_right = None

        if params.output_file_left:
            output_file_left = open(params.output_file_left, 'w')

        if params.output_file_right:
            output_file_right = open(params.output_file_right, 'w')

        if params.output_stats_file_left:
            output_stats_file_left = open(params.output_stats_file_left, 'w')

        if params.output_stats_file_right:
            output_stats_file_right = open(params.output_stats_file_right, 'w')

        mi = ['L']
        if params.diploid:
            mi = ['L', 'R']

        logger.warning(f'Processing Chromosome {params.chromosome}...')

        iterators = []
        discard_functions = []

        tabix = False

        if tabix:

            for i, file_info in enumerate(params.vcf_files):
                try:
                    vcf_tabix = pysam.TabixFile(file_info.file_name)
                except OSError:
                    # attempt CSI index if TBI can't be found
                    vcf_tabix = pysam.TabixFile(
                        filename=file_info.file_name,
                        index=f'{file_info.file_name}.csi',
                    )

                try:
                    vcf_iterator = vcf_tabix.fetch(
                        params.chromosome,  # 3052782, 3052784,
                        parser=pysam.asVCF(),
                    )
                    iterators.append(vcf_iterator)
                except ValueError:
                    iterators.append(None)

                if file_info.discard_file:
                    vcf_discard = open(file_info.discard_file, 'w')

                    def discard_record(rec):
                        vcf_discard.write(str(rec))
                        vcf_discard.write('\n')

                    discard_functions.append(discard_record)
                else:
                    discard_functions.append(lambda rec: None)

        else:
            for i, file_info in enumerate(params.vcf_files):
                try:
                    vcf_vf = pysam.VariantFile(file_info.file_name)
                except OSError:
                    # attempt CSI index if TBI can't be found
                    vcf_vf = pysam.VariantFile(
                        filename=file_info.file_name,
                        index_filename=f'{file_info.file_name}.csi',
                    )

                try:
                    vcf_iterator = vcf_vf.fetch(
                        params.chromosome  # , 3052780, 3052785
                    )
                    iterators.append(vcf_iterator)
                except ValueError:
                    iterators.append(None)

                if file_info.discard_file:
                    vcf_discard = open(file_info.discard_file, 'w')

                    def discard_record(rec):
                        vcf_discard.write(str(rec))
                        vcf_discard.write('\n')

                    discard_functions.append(discard_record)
                else:
                    discard_functions.append(lambda rec: None)

        n = 0

        transcript_info = None
        transcript_features = None
        transcript_tree = None
        transcript_info_left = None
        transcript_info_right = None

        if params.gtf_file:
            transcript_info = gtf.parse_gtf(params.gtf_file, params.chromosome)
            transcript_features = transcript_info['features']
            transcript_tree = transcript_info['tree']

            transcript_info_left = copy.deepcopy(transcript_features)
            transcript_info_right = copy.deepcopy(transcript_features)

        line_numbers = 0
        # print('iterators=' + str(type(iterators)))
        # print('iterators[0]=' + str(type(iterators[0])))
        for vcf_records in walk_vcfs_together(iterators):
            # print('vcf_records=' + str(type(vcf_records)))
            for i, vcf_record in enumerate(vcf_records):
                # print('vcf_record=' + str(type(vcf_record)))
                # logger.debug(vcf_record)
                if vcf_record is None:
                    continue
                # logger.debug(vcf_record.alt)
                # logger.debug(type(vcf_record.alt))
                logger.debug('------------')
                # print(f'{vcf_record.pos=}')

                if tabix:
                    gt = vcf.parse_gt_tuple_orig(
                        vcf_record, params.vcf_files[i].sample_index
                    )
                else:
                    gt = vcf.parse_gt_tuple(
                        vcf_record, params.vcf_files[i].sample_name
                    )

                # logger.debug(gt)
                logger.debug(f'{gt=}')

                line_numbers = line_numbers + 1
                if gt.is_snp:
                    # snp
                    if params.passed and 'PASS' not in vcf_record.filter:
                        discard_functions[i](vcf_record)

                    # logger.debug(f"Processing SNP {vcf_record}")
                    n += 1

                    if params.quality and gt.fi == '0':
                        discard_functions[i](vcf_record)
                    elif gt.left is None or gt.right is None:
                        discard_functions[i](vcf_record)
                    else:
                        if params.diploid:
                            # 0 is the same as REF and do not need
                            if gt.gt_left != 0:
                                vpos = vcf_record.pos + (1 if tabix else 0)

                                output_file_left.write(
                                    f'{params.chromosome}_L\t'
                                    f'{vpos}\t'
                                    '.\t'
                                    f'{vcf_record.ref}\t'
                                    f'{gt.left}\t'
                                    '.\n'
                                )

                                if params.gtf_file:
                                    overlapping_features = transcript_tree[
                                        vpos
                                    ]
                                    for interval in overlapping_features:
                                        feature_id = interval.data
                                        transcript_info_left[
                                            feature_id
                                        ].snp_count += 1

                            if gt.gt_right != 0:
                                vpos = vcf_record.pos + (1 if tabix else 0)

                                output_file_right.write(
                                    f'{params.chromosome}_R\t'
                                    f'{vpos}\t'
                                    '.\t'
                                    f'{vcf_record.ref}\t'
                                    f'{gt.right}\t'
                                    '.\n'
                                )

                                if params.gtf_file:
                                    overlapping_features = transcript_tree[
                                        vpos
                                    ]
                                    for interval in overlapping_features:
                                        feature_id = interval.data
                                        transcript_info_right[
                                            feature_id
                                        ].snp_count += 1

                        else:
                            if gt.gt_left == gt.gt_right and gt.gt_left != 0:
                                # ignore heterozygotes 0/1, 1/0,
                                #     only process 0/0 and 1/1
                                # logger.debug('ACCEPTED')
                                # logger.debug(
                                #     f'pos {vcf_snp.pos} : '
                                #     f'ref {vcf_snp.ref}, '
                                #     f'left {gt.left},
                                #     f'right {gt.right}'
                                # )
                                vpos = vcf_record.pos + (1 if tabix else 0)

                                output_file_left.write(
                                    f'{params.chromosome}\t'
                                    f'{vpos}\t'
                                    '.\t'
                                    f'{vcf_record.ref}\t'
                                    f'{gt.left}\t'
                                    '.\n'
                                )

                                if params.gtf_file:
                                    overlapping_features = transcript_tree[
                                        vpos
                                    ]
                                    for interval in overlapping_features:
                                        feature_id = interval.data
                                        transcript_info_left[
                                            feature_id
                                        ].snp_count += 1
                else:
                    # indel
                    logger.debug(f'Processing INDEL {vcf_record}')

                    if params.passed and 'PASS' not in vcf_record.filter:
                        logger.debug('TOSSED: FILTERED ON PASS')
                        logger.debug(vcf_record)
                        stats = update_stats(stats, 'FILTERED ON PASS')
                        discard_functions[i](vcf_record)
                        continue

                    elif params.quality and gt.fi == '0':
                        # FI : Whether a sample was a Pass(1) or
                        #      fail (0) based on FILTER values

                        logger.debug('TOSSED: FILTERED ON QUALITY')
                        logger.debug(vcf_record)
                        stats = update_stats(stats, 'FILTERED ON QUALITY')
                        discard_functions[i](vcf_record)
                        continue

                    elif gt.left is None and gt.right is None:
                        logger.debug('TOSSED: NO STRAIN DATA')
                        logger.debug(vcf_record)
                        stats = update_stats(stats, 'NO STRAIN DATA')
                        # logger.debug(i)
                        # logger.debug(type(vcf_record))
                        discard_functions[i](vcf_record)
                        continue

                    elif not params.diploid and gt.left != gt.right:
                        # haploid or hexaploid
                        # gt must be equal
                        logger.debug('TOSSED: HETEROZYGOUS')
                        logger.debug(vcf_record)
                        stats = update_stats(stats, 'HETEROZYGOUS')
                        discard_functions[i](vcf_record)
                        continue

                    # START L AND R, ONLY R IF DIPLOID

                    for l_or_r in mi:
                        # logger.debug('******************')
                        # logger.debug(l_or_r)
                        if l_or_r == 'L':
                            # logger.debug('->LEFT')
                            lr_out = '_L' if params.diploid else ''
                            alt_seq = str(gt.left)
                            stats = params.stats_left
                            output_file = output_file_left
                            prev_next_ref_pos = params.prev_next_ref_pos_left
                        else:
                            # logger.debug('->RIGHT')
                            lr_out = '_R' if params.diploid else ''
                            alt_seq = str(gt.right)
                            stats = params.stats_right
                            output_file = output_file_right
                            prev_next_ref_pos = params.prev_next_ref_pos_right

                        logger.debug(f'prev_next_ref_pos={prev_next_ref_pos}')

                        if gt.ref == alt_seq:
                            logger.debug('TOSSED, REF AND ALT ARE EQUAL')
                            logger.debug(vcf_record)
                            stats = update_stats(
                                stats, 'REF AND ALT ARE EQUAL'
                            )
                            discard_functions[i](vcf_record)
                            continue

                        orig_alt_seq = alt_seq

                        # s = vcf_record[
                        #    params.vcf_files[i].sample_index
                        # ]
                        # logger.debug(f'SAMPLE: {s}')
                        logger.debug(
                            f'REF="{gt.ref}", ALT_L="{gt.left}", '
                            f'ALT_R="{gt.right}", POS={vcf_record.pos}'
                        )

                        position = vcf_record.pos + (1 if tabix else 0)

                        ref_seq = str(gt.ref)
                        len_ref = len(ref_seq)
                        len_alt = len(alt_seq)

                        base_pos_diff = 0

                        if position < prev_next_ref_pos:
                            logger.debug(f'TOSSED: VCF ROLLBACK: {vcf_record}')
                            logger.debug(vcf_record)

                            stats = update_stats(stats, 'VCF ROLLBACK')
                            discard_functions[i](vcf_record)
                            continue

                        # find the position where the first base change is
                        for n in range(min(len_ref, len_alt)):
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

                        # if dt > dq:
                        #    logger.debug("         DELETION:")
                        #    next_ref_pos = position
                        # else:
                        #    logger.debug("        INSERTION:")
                        #    next_ref_pos = position + len(ref_seq)

                        fragment_size = position - prev_next_ref_pos
                        base_changes = dq - dt

                        logger.debug(f'           gt.ref: {gt.ref}')
                        logger.debug(f'          ref_seq: {ref_seq}')
                        logger.debug(f'               dt: {dt}')
                        logger.debug(f'           gt.alt: {orig_alt_seq}')
                        logger.debug(f'          alt_seq: {alt_seq}')
                        logger.debug(f'               dq: {dq}')
                        logger.debug(f'         position: {position}')
                        logger.debug(f'prev_next_ref_pos: {prev_next_ref_pos}')
                        logger.debug(f'     next_ref_pos: {next_ref_pos}')
                        logger.debug(f'    fragment_size: {fragment_size}')
                        logger.debug(f'     base_changes: {base_changes}')
                        logger.debug(f'    base_pos_diff: {base_pos_diff}')
                        logger.debug(f'     shared_bases: {shared_bases}')

                        # fix any 0 length
                        if fragment_size < 0:
                            # logger.debug(f'TOSSED: FRAGMENT: {vcf_record}')

                            stats = update_stats(stats, 'FRAGMENT SIZE < 0')
                            discard_functions[i](vcf_record)
                            continue

                        if fragment_size != 0:
                            ref_str = ref_seq if ref_seq else '.'
                            alt_str = alt_seq if alt_seq else '.'
                            out = (
                                f'{params.chromosome}{lr_out}\t'
                                f'{vcf_record.pos + (1 if tabix else 0)}\t'
                                f'{shared_bases}\t{ref_str}\t{alt_str}\t'
                                f'{fragment_size}\n'
                            )
                            logger.debug(out)
                            output_file.write(out)

                            vpos = vcf_record.pos + (1 if tabix else 0)
                            overlapping_features = transcript_tree[vpos]
                            if l_or_r == 'L':
                                for interval in overlapping_features:
                                    feature_id = interval.data
                                    transcript_info_left[
                                        feature_id
                                    ].indel_count += 1
                            else:
                                for interval in overlapping_features:
                                    feature_id = interval.data
                                    transcript_info_right[
                                        feature_id
                                    ].indel_count += 1

                        else:
                            # THIS SHOULD NOT HAPPEN
                            raise G2GVCFError('Conflicting VCF entries')

                        stats = update_stats(stats, 'ACCEPTED')

                        if l_or_r == 'L':
                            params.stats_left = stats
                            params.prev_next_ref_pos_left = next_ref_pos
                            logger.debug(
                                'setting params.prev_next_ref_pos_left='
                                f'{params.prev_next_ref_pos_left}'
                            )
                        else:
                            params.stats_right = stats
                            params.prev_next_ref_pos_right = next_ref_pos
                            logger.debug(
                                'setting params.prev_next_ref_pos_right='
                                f'{params.prev_next_ref_pos_right}'
                            )

        if params.output_file_left:
            output_file_left.close()

        if params.output_file_right:
            output_file_right.close()

        # output all gene info from stats
        if params.output_stats_file_left:
            transcript_list = list(transcript_info_left.values())

            # Sort the list by gene_id first, then by transcript_id
            sorted_transcripts = sorted(
                transcript_list, key=lambda x: (x.gene_id, x.transcript_id)
            )

            for transcript in sorted_transcripts:
                output_stats_file_left.write(
                    f'{transcript.gene_id}\t'
                    f'{transcript.gene_name}\t'
                    f'{transcript.transcript_id}\t'
                    f'{transcript.snp_count}\t'
                    f'{transcript.indel_count}\n'
                )

            output_file_left.close()

        if params.output_stats_file_right:
            transcript_list = list(transcript_info_right.values())

            # Sort the list by gene_id first, then by transcript_id
            sorted_transcripts = sorted(
                transcript_list, key=lambda x: (x.gene_id, x.transcript_id)
            )

            for transcript in sorted_transcripts:
                output_stats_file_right.write(
                    f'{transcript.gene_id}\t'
                    f'{transcript.gene_name}\t'
                    f'{transcript.transcript_id}\t'
                    f'{transcript.snp_count}\t'
                    f'{transcript.indel_count}\n'
                )

            output_file_right.close()

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        g2g_utils.show_error()
        raise Exception(f'Unknown exception: {e}')

    return {
        'chrom': params.chromosome,
        'stats': stats,
        'params': params,
        'line_numbers': line_numbers,
    }


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    Args:
        args: The arguments to process_piece.

    Returns:
        The same as process_piece
    """
    try:
        result=process_piece(*args)
    except Exception as e:
        print (args)
        raise e
    return result


def create_vci_header(
    temp_directory: str,
    fasta_file: str,
    vcf_input_files: list[VCFFileInformation],
    vci_contigs: list[str],
    strain: str,
    vcf_keep: bool,
    passed: bool,
    quality: bool,
    diploid: bool,
    num_processes: int,
) -> str:
    """
    Create a VCI file header that contains meta information about the VCI file.

    Args:
        temp_directory: Directory to build the file.
        fasta_file: Fasta file used in VCI file creation.
        vcf_input_files: List of VCIFileInformation.
        vci_contigs: Ordered list of contigs.
        strain: The strain used to make the file.
        vcf_keep: True to place troubling VCF lines in extra file.
        passed: True uses only VCF lines that have a PASS for the filter.
        quality: True to filter on quality, FI=PASS.
        diploid: Create diploid VCI file.
        num_processes: Specify the number of processes.

    """
    file = g2g_utils.gen_file_name(
        'header', output_dir=temp_directory, extension='', append_time=False
    )

    with open(file, 'w') as fd:
        create_time = time.strftime('%m/%d/%Y %H:%M:%S')
        fd.write(f'##CREATION_TIME={create_time}\n')

        for vcf_file in vcf_input_files:
            fd.write(f'##INPUT_VCF={vcf_file.file_name}\n')

        fd.write(f'##FASTA_FILE={fasta_file}\n')
        fd.write(f'##STRAIN={strain}\n')
        fd.write(f'##VCF_KEEP={vcf_keep}\n')
        fd.write(f'##FILTER_PASSED={passed}\n')
        fd.write(f'##FILTER_QUALITY={quality}\n')
        fd.write(f'##DIPLOID={diploid}\n')
        fd.write(f'##PROCESSES={num_processes}\n')

        fasta_file = fasta.FastaFile(fasta_file)

        for contig in vci_contigs:
            fd.write(
                f'##CONTIG={contig}:{fasta_file.get_reference_length(contig)}\n'
            )

        fd.write('#CHROM\tPOS\tANCHOR\tDEL\tINS\tFRAG\n')
    return file


def process(
    vcf_files: list[str],
    fasta_file: str,
    output_file: str,
    strain: str,
    gtf_file: str = None,
    vcf_keep: bool = False,
    passed: bool = False,
    quality: bool = False,
    diploid: bool = False,
    num_processes: int | None = None,
    bgzip: bool = False,
) -> None:
    """
    Parse the VCF file and create a VCI file.

    Args
        vcf_files: Name of the VCF files.
        fasta_file: Name of the Fasta file associated with the VCF file.
        output_file: Name of output file.
        strain: Which strain to process.
        gtf_file: Optional GTF file for additional annotations.
        vcf_keep: True to place troubling VCF lines in extra file.
        passed: True uses only VCF lines that have a PASS for the filter.
        quality: True to filter on quality, FI=PASS.
        diploid: Create diploid VCI file.
        num_processes: Specify the number of processes.
        bgzip: True tp bgzip the VCI file.
    """
    start = time.time()

    output_file = g2g_utils.check_file(output_file, 'w')
    output_file_dir = os.path.dirname(output_file)
    fasta_file = g2g_utils.check_file(fasta_file)

    output_stats_file = None

    if gtf_file:
        gtf_file = g2g_utils.check_file(gtf_file)
        output_stats_file = f'{output_file}.stats.tsv'

    vcf_file_inputs = []

    if vcf_files:
        for file_name in vcf_files:
            vcf_file = g2g_utils.check_file(file_name)
            logger.warning(f'Input VCF file: {vcf_file}')
            logger.debug('Checking for index file, creating if needed...')
            g2g_utils.index_file(
                file_name=vcf_file, file_format='vcf', overwrite=False
            )

            vcf_discard_file = None
            if vcf_keep:
                vcf_discard_file = os.path.join(
                    output_file_dir,
                    f'{os.path.basename(output_file)}.errors.vcf',
                )
                logger.warning(
                    f'Output VCF indel discard file: {vcf_discard_file}'
                )

            vcf_file_inputs.append(
                VCFFileInformation(vcf_file, vcf_discard_file)
            )

    if len(vcf_file_inputs) == 0:
        raise G2GValueError('No VCF files.')

    if not fasta_file:
        raise G2GValueError('No fasta file was specified.')

    if not strain:
        raise G2GValueError('No strain was specified.')

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    logger.warning(f'Input Fasta File: {fasta_file}')
    logger.warning(f'Strain: {strain}')
    logger.warning(f'Pass filter on: {passed}')
    logger.warning(f'Quality filter on: {quality}')
    logger.warning(f'Diploid: {diploid}')
    logger.warning(f'BGZip: {bgzip}')
    logger.debug(f'Number of processes: {num_processes}')

    if gtf_file:
        logger.warning(f'GTF File: {gtf_file}')
        logger.warning(f'Stats File: {output_stats_file}')

    if bgzip:
        logger.warning(f'Output VCI File: {output_file}.gz')
    else:
        logger.warning(f'Output VCI File: {output_file}')

    # not all chromosomes/seq_id will be processed if not in vcf file
    processed_seq_ids = {}

    temp_directory = g2g_utils.create_temp_dir('vcf2vci', dir='.')
    logger.debug(f'Temp directory: {temp_directory}')

    for i, vcf_file in enumerate(vcf_file_inputs):
        try:
            tb_file = pysam.TabixFile(vcf_file.file_name)
        except OSError:
            # attempt CSI index if TBI can't be found
            tb_file = pysam.TabixFile(
                filename=vcf_file.file_name, index=f'{vcf_file.file_name}.csi'
            )
        for h in tb_file.header:
            h = g2g_utils.s(h)
            if h[:6] == '#CHROM':
                try:
                    elem = h.split('\t')
                    samples = elem[9:]
                    samples = dict(
                        zip(samples, (x for x in range(len(samples))))
                    )

                    vcf_file_inputs[i].sample_name = strain
                    vcf_file_inputs[i].sample_index = samples[strain]
                except KeyError:
                    elem = h.split('\t')
                    valid_strains = ', '.join(elem[9:])
                    raise G2GVCFError(
                        f'Unknown strain "{strain}", '
                        f'valid strains are: {valid_strains}'
                    )

        for seq_id in tb_file.contigs:
            processed_seq_ids[seq_id] = False

    tmp_processed_seq_ids = OrderedDict()
    vci_file_contigs = []

    for k in g2g_utils.natsorted(processed_seq_ids.keys()):
        tmp_processed_seq_ids[k] = False
        vci_file_contigs.append(k)

    processed_seq_ids = tmp_processed_seq_ids

    header_file = create_vci_header(
        temp_directory,
        fasta_file,
        vcf_file_inputs,
        vci_file_contigs,
        strain,
        vcf_keep,
        passed,
        quality,
        diploid,
        num_processes,
    )

    all_params = []
    pool = None

    try:
        for c in processed_seq_ids:
            params = VCF2VCIInfo()
            params.chromosome = c
            params.vcf_files = vcf_file_inputs
            params.fasta_file = fasta_file
            params.gtf_file = gtf_file
            params.diploid = diploid
            params.passed = passed
            params.quality = quality
            params.vcf_keep = vcf_keep
            params.log_level = logger.level

            if diploid:
                params.output_file_left = g2g_utils.gen_file_name(
                    f'chr{c}.left',
                    output_dir=temp_directory,
                    extension='vci',
                    append_time=False,
                )

                params.output_file_right = g2g_utils.gen_file_name(
                    f'chr{c}.right',
                    output_dir=temp_directory,
                    extension='vci',
                    append_time=False,
                )

                # delete in case the files already exist from a previous run
                g2g_utils.delete_file(params.output_file_left)
                g2g_utils.delete_file(params.output_file_right)

                if gtf_file:
                    params.output_stats_file_left = g2g_utils.gen_file_name(
                        f'chr{c}.left.stats',
                        output_dir=temp_directory,
                        extension='tsv',
                        append_time=False,
                    )

                    params.output_stats_file_right = g2g_utils.gen_file_name(
                        f'chr{c}.right.stats',
                        output_dir=temp_directory,
                        extension='tsv',
                        append_time=False,
                    )

                    g2g_utils.delete_file(params.output_stats_file_left)
                    g2g_utils.delete_file(params.output_stats_file_right)

            else:
                params.output_file_left = g2g_utils.gen_file_name(
                    f'chr{c}.left',
                    output_dir=temp_directory,
                    extension='vci',
                    append_time=False,
                )

                # delete in case the files already exist from a previous run
                g2g_utils.delete_file(params.output_file_left)

                if gtf_file:
                    params.output_stats_file_left = g2g_utils.gen_file_name(
                        f'chr{c}.left.stats',
                        output_dir=temp_directory,
                        extension='tsv',
                        append_time=False,
                    )

                    g2g_utils.delete_file(params.output_stats_file_left)

            all_params.append(params)

        logger.info('Parsing VCF files...')

        args = zip(all_params)
        with multiprocessing.Pool(num_processes) as pool:
            results = pool.map(wrapper, args)

        # parse results
        total = 0
        for r in results:
            total += r['line_numbers']

        logger.debug('Combining temp files...')
        logger.warning('Finalizing VCI File')
        logger.warning(f'Parsed {total:,} lines')

        files = [header_file]

        for mi in all_params:
            files.append(mi.output_file_left)
            if diploid:
                files.append(mi.output_file_right)

        g2g_utils.concatenate_files(files, output_file, True)

        if gtf_file:
            logger.warning('Finalizing Stats File')
            files_stats_left = []
            files_stats_right = []

            if diploid:
                output_stats_file_left = g2g_utils.adjust_file_name(
                    output_stats_file, 'L'
                )
                output_stats_file_right = g2g_utils.adjust_file_name(
                    output_stats_file, 'R'
                )
            else:
                output_stats_file_left = output_stats_file
                output_stats_file_right = None

            for mi in all_params:
                files_stats_left.append(mi.output_stats_file_left)
                if diploid:
                    files_stats_right.append(mi.output_stats_file_right)

            g2g_utils.concatenate_files(
                files_stats_left, output_stats_file_left, True
            )

            if diploid:
                g2g_utils.concatenate_files(
                    files_stats_right, output_stats_file_right, True
                )

        if bgzip:
            logger.warning('Compressing VCI File and indexing')
            g2g_utils.bgzip_and_index_file(
                output_file,
                f'{output_file}.gz',
                delete_original=True,
                file_format='vcf',
            )

        logger.warning('VCI File created')

    except KeyboardInterruptError:
        pool.terminate()
        raise G2GError('Keyboard quit consumed')
    except KeyboardInterrupt:
        pool.terminate()
        raise G2GError('Execution halted')
    except Exception as e:
        g2g_utils.show_error()
        pool.terminate()
        raise G2GError(f'Execution halted unknown error: {e}')
    finally:
        g2g_utils.delete_dir(temp_directory)
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')

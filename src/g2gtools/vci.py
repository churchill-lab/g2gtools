# standard library imports
import collections
import time
from collections import OrderedDict

# 3rd party library imports
import pysam
from bx.intervals.intersection import Interval, IntervalTree

# local library imports
from g2gtools.exceptions import G2GError
from g2gtools.exceptions import G2GRegionError
import g2gtools.g2g_utils as g2g_utils
import g2gtools.region as region

logger = g2g_utils.get_logger('g2gtools')


InfoFields = ['chr', 'start', 'end', 'shared', 'inserted', 'deleted', 'pos']
IntervalInfo = collections.namedtuple('IntervalInfo', InfoFields)

MappingFields = [
    'from_chr',
    'from_start',
    'from_end',
    'from_seq',
    'to_chr',
    'to_start',
    'to_end',
    'to_seq',
    'same_bases',
    'vcf_pos',
]
IntervalMapping = collections.namedtuple('IntervalMapping', MappingFields)

#
# VCIFile
#
# HEADERS
# ##EXECUTION_TIME
# ##INDEL_VCF
# ##INDEL_SNP
# ##STRAIN
# ##VCF_KEEP
# ##FILTER_PASSED
# ##FILTER_QUALITY
# ##DIPLOID
# ##PROCESSES
# #CHROM POS SHARE INS DEL FRAG
# DATA


def intersect_regions(
    chr1: str, start1: int, end1: int, chr2: str, start2: int, end2: int
) -> tuple[str, int, int] | None:
    """
    Return the intersection of two regions.

    Args:
        chr1: First chromosome.
        start1: First start position.
        end1: First end position.
        chr2: Second chromosome.
        start2: Second start position.
        end2: Second end position.

    Returns:
        The intersecting region as a tuple of chromosome, start, end or
        None if the regions do not overlap.

    Raises:
        G2GRegionError: If the region is invalid.
    """

    if chr1 != chr2:
        return None

    if int(start1) > int(end2) or int(end1) < int(start2):
        return None

    if int(start1) > int(end1) or int(start2) > int(end2):
        raise G2GRegionError('Start cannot be larger than end')

    return chr1, max(start1, start2), min(end1, end2)


class VCIFile:
    """
    Encapsulate a VCI (Variant Call Information) File
    """

    HEADERS = [
        'CREATION_TIME',
        'INDEL_VCF',
        'INDEL_SNP',
        'STRAIN',
        'VCF_KEEP',
        'FILTER_PASSED',
        'FILTER_QUALITY',
        'DIPLOID',
        'PROCESSES',
    ]

    def __init__(
        self,
        file_name: str,
        mode: str = 'r',
        parser: pysam.libctabix.Parser | None = None,
        index: str | None = None,
        encoding: str = 'ascii',
        seq_ids: list[str] | None = None,
        reverse: bool = False,
    ):
        """
        Initialize VCI File.

        Args:
            file_name: The VCI file name.
            mode: Mode to read the file.
            parser: Sets the default parser. If None, the results are returned
                as an un-parsed string.
            index: The filename of the index. If not set, the default is to
                assume that the index is called.
            encoding: The encoding passed to the parser.
            seq_ids: Override what is in the file and just use these seq ids.
            reverse: Is the file reversed.
        """
        self.file_name: str = file_name
        try:
            self.file_name = self.file_name.decode()
        except (UnicodeDecodeError, AttributeError):
            pass

        self.mapping_tree: dict[str, IntervalTree] = {}
        self.headers: dict[str, str] = {}
        self.seq_ids: list[str] = seq_ids
        self.is_reversed: bool = reverse
        self.contigs: dict[str, int] = {}
        self.valid: bool = False
        self.debug_level: int = 0

        self.tabix_file: pysam.TabixFile = pysam.TabixFile(
            self.file_name,
            mode=mode,
            parser=parser,
            index=index,
            encoding=encoding,
        )

        g2g_utils.index_file(
            file_name=file_name, file_format='vci', overwrite=False
        )

        self.parse_header()

    def get_filename(self) -> str:
        """
        Get the filename of the VCI File.

        Returns:
            The filename.
        """
        return self.file_name

    def __getattr__(self, name):
        """
        Get the attribute of the VCI File.

        Args:
            name: The attribute name.

        Returns:
            The attribute value.
        """
        return getattr(self.tabix_file, name)

    def fetch(
        self,
        reference: str | None = None,
        start: int | None = None,
        end: int | None = None,
        region: str | None = None,
        parser: pysam.libctabix.Parser | None = None,
        multiple_iterators: bool | None = False,
    ) -> pysam.libctabix.TabixIterator:
        """
        Fetch reads in the specified region or reference. Without a contig or
        region all mapped reads in the file will be fetched.

        Args:
            reference: The chromosome or contig.
            start: Start position.
            end: End position.
            region: A region string in samtools format.
            parser: Parser to use.
            multiple_iterators: Whether multiple iterators can be used.

        Returns:
            An iterator over a collection of reads.
        """
        # todo: potentially alter the reference value or the region value???
        return self.tabix_file.fetch(
            reference, start, end, region, parser, multiple_iterators
        )

    def is_diploid(self):
        """
        Check if this is a DIPLOID VCI File.

        Returns:
            True if this is a DIPLOID VCI File.
        """
        if 'DIPLOID' in self.headers:
            return self.headers['DIPLOID'].lower() in ['true', 't', 'yes', '1']
        return False

    def is_haploid(self):
        """
        Check if this is a HAPLOID VCI File.

        Returns:
            True if this is a HAPLOID VCI File.
        """
        if 'DIPLOID' in self.headers:
            if self.headers['DIPLOID'].lower() in ['true', 't', 'yes', '1']:
                return False
        return True

    def parse_header(self):
        """
        Parse the VCI header file.
        """
        self.headers = {}

        for header_line in self.tabix_file.header:
            header_line = g2g_utils.s(header_line)

            if header_line.startswith('##'):
                elem = header_line[2:].split('=')

                if elem[0] == 'CONTIG':
                    info = elem[1].split(':')
                    self.contigs[info[0]] = int(info[1])
                else:
                    self.headers[elem[0]] = elem[1]

    def get_seq_ids(self):
        """
        Return all the contigs or chromosomes in the mapping tree.

        Returns:
            A list of the keys of the mapping tree.
        """
        if self.mapping_tree:
            return list(self.mapping_tree.keys())
        return None

    def parse(self, reverse: bool) -> None:
        """
        Parse the VCI File.

        Args:
            reverse: Whether to switch inserted and deleted.
        """
        # start = time.time()
        mapping_tree = {}
        try:
            total_num_lines_chrom = 0
            total_num_lines_processed = 0

            _chrom = 0
            _pos = 1
            _shared = 2
            _deleted = 3 if not reverse else 4
            _inserted = 4 if not reverse else 3
            _fragment = 5

            contigs = self.tabix_file.contigs

            if self.seq_ids:
                contigs = self.seq_ids

            # create a bx tree for the indels
            for contig in contigs:
                # contig_start_time = time.time()

                num_lines_chrom = 0
                num_lines_processed = 0

                pos_from = 0
                pos_to = 0

                if contig not in mapping_tree:
                    mapping_tree[contig] = IntervalTree()

                iterator = None

                try:
                    iterator = self.tabix_file.fetch(
                        contig, parser=pysam.asTuple()
                    )
                except Exception:
                    pass

                if iterator is None:
                    continue

                for rec in iterator:
                    num_lines_chrom += 1
                    total_num_lines_chrom += 1

                    if len(rec) != 6:
                        raise G2GError(
                            'Unexpected line in VCI file. Line #'
                            f'{total_num_lines_chrom:,}: {rec}'
                        )

                    if rec[2] == '.':
                        continue

                    if rec[_inserted] == '.':
                        inserted_bases = 0
                    else:
                        inserted_bases = len(rec[_inserted])

                    if rec[_deleted] == '.':
                        deleted_bases = 0
                    else:
                        deleted_bases = len(rec[_deleted])

                    fragment = int(rec[_fragment])
#g2gtools [09:35:22] Inserting interval 195349534 - 195359214
#g2gtools [09:35:22] pos_from=195359216, pos_to=195305722
#g2gtools [09:35:22] Inserting interval 195359216 - 195363317
# 1	195363316	CA	A	.	4101
# 1	195365691	.	G	T	.


                    # logger.debug(f'pos_from={pos_from}, pos_to={pos_to}')
                    # logger.debug(
                    #      f'Inserting interval {pos_from} - {pos_from+fragment}'
                    # )
                    # 'chr', 'start', 'end', 'shared', 'inserted', 'deleted', 'pos'
                    info = IntervalInfo(
                        contig,  # chr
                        pos_to,  # start
                        pos_to + fragment,  # end
                        rec[_shared],  # shared
                        rec[_inserted],  # inserted
                        rec[_deleted],  # deleted
                        rec[_pos],  # pos
                    )
                    interval = Interval(pos_from, pos_from + fragment, info)

                    # logging.debug(interval)

                    mapping_tree[contig].insert_interval(interval)

                    pos_from += fragment + deleted_bases
                    pos_to += fragment + inserted_bases
                    num_lines_processed += 1
                    total_num_lines_processed += 1

                if self.is_diploid():
                    info = IntervalInfo(
                        contig,  # chr
                        pos_to,  # start
                        pos_to + (self.contigs[contig[:-2]] - pos_from),  # end
                        None,  # shared
                        None,  # inserted
                        None,  # deleted
                        None,  # pos
                    )
                    interval = Interval(
                        pos_from, self.contigs[contig[:-2]], info
                    )
                else:
                    info = IntervalInfo(
                        contig,  # chr
                        pos_to,  # start
                        pos_to + (self.contigs[contig] - pos_from),  # end
                        None,  # shared
                        None,  # inserted
                        None,  # deleted
                        None,  # pos
                    )
                    interval = Interval(pos_from, self.contigs[contig], info)

                mapping_tree[contig].insert_interval(interval)

                # elapsed_time = g2g_utils.format_time(
                #     contig_start_time, time.time()
                # )
                # logging.debug(
                #     f'Parsed {num_lines_processed:,} '
                #     f'lines for contig {contig} in {elapsed_time}'
                # )

            self.valid = True
            self.mapping_tree = mapping_tree
        except Exception:
            g2g_utils.show_error()

        # fmt_time = g2g_utils.format_time(start, time.time())
        # logging.info('Parsing complete: {fmt_time}')

    def add_to_tree(self, contig: str, tree: IntervalTree) -> None:
        """
        Add tree to the mapping tree.

        Args:
            contig: The contig.
            tree: The tree to add.
        """
        if contig not in self.mapping_tree:
            self.mapping_tree[contig] = tree

    def find_mappings(
        self, chromosome: str, start: int, end: int
    ) -> list[IntervalMapping] | None:
        """
        Find mapping from source to target.

        Args:
            chromosome: The chromosome.
            start: The start position.
            end: The end position.

        Returns:
            A list on IntervalMappings or None if none are found.
        """
        mappings = []

        if chromosome not in self.mapping_tree:
            # logging.debug(f'Chromosome {chromosome} not in mapping tree')
            # logging.debug(f'Chromosomes: {list(self.mapping_tree.keys())}')
            return None

        all_intervals = self.mapping_tree[chromosome].find(start, end)

        if len(all_intervals) == 0:
            # logging.debug('No intervals found')
            return None
        else:
            for interval in all_intervals:
                # logging.debug(f'Interval {len(mappings)} = {interval}')
                chromosome, real_start, real_end = intersect_regions(
                    chromosome,
                    start,
                    end,
                    chromosome,
                    interval.start,
                    interval.end,
                )
                offset = abs(real_start - interval.start)
                size = abs(real_end - real_start)
                i_start = interval.value[1] + offset

                # 'from_chr', 'from_start', 'from_end', 'from_seq',
                # 'to_chr', 'to_start', 'to_end', 'to_seq',
                # 'same_bases', 'vcf_pos'

                mappings.append(
                    IntervalMapping(
                        chromosome,
                        real_start,
                        real_end,
                        interval.value.inserted,
                        interval.value.chr,
                        i_start,
                        i_start + size,
                        interval.value.deleted,
                        interval.value.shared,
                        interval.value.pos,
                    )
                )
                # logging.debug(f'Mapping {len(mappings)-1}={mappings[-1]}')

        return mappings


def vci_query(vci_file_name_in: str, reg: region.Region) -> None:
    """
    Query the VCI file and output the matching records.

    Args:
        vci_file_name_in: The name of the VCI file.
        region: The region to query.
    """
    start = time.time()

    vci_file_name_in = g2g_utils.check_file(vci_file_name_in, 'r')

    logger.warning(f'VCI File: {vci_file_name_in}')
    logger.warning(f'Region: {reg}')
    logger.debug(f'seq_id: {reg.seq_id}')
    logger.debug(f'start: {reg.start}')
    logger.debug(f'end: {reg.end}')

    vci_f = VCIFile(vci_file_name_in, seq_ids=[reg.seq_id])
    vci_f.parse(False)

    mappings = vci_f.find_mappings(reg.seq_id, reg.start, reg.end)

    if mappings is None:
        logger.warning('No mappings found')
        return

    for m in mappings:
        logger.debug(m)

    start_pos = mappings[0].to_start
    end_pos = mappings[-1].to_end

    logger.debug(f'Converted: {reg.seq_id}:{start_pos + 1}-{end_pos + 1}')

    for line in vci_f.fetch(
        reference=reg.seq_id,
        start=start_pos,
        end=end_pos,
        parser=pysam.asTuple(),
    ):
        print(str(line))

    fmt_time = g2g_utils.format_time(start, time.time())
    logger.warning(f'VCI Query Complete: {fmt_time}')


def convert_region(
        vci_file: str | VCIFile,
        reg: region.Region | list[region.Region] | str | list[str]
    ):
    """
    Convert location(s) to new coordinates.

    :param vci_file: a string or VCIFile
    :param reg: a location or list of locations
    :return: a list of regions or a single region
    """

    all_regions = []

    if isinstance(reg, region.Region):
        all_regions.append(region)
    elif isinstance(reg, str):
        all_regions.append(region.parse_region(region))
    elif isinstance(reg, list):
        for loc in reg:
            if isinstance(loc, region.Region):
                all_regions.append(loc)
            elif isinstance(loc, str):
                all_regions.append(region.parse_region(loc))
            else:
                raise G2GRegionError(f'Unknown region type: {type(loc)}')
    else:
        raise G2GRegionError(f'Unknown region type: {type(region)}')

    if isinstance(vci_file, VCIFile):
        logger.warning(f'VCI FILE: {vci_file.filename}')
        logger.info(f'VCI FILE IS DIPLOID: {vci_file.is_diploid()}')
    else:
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = VCIFile(vci_file)
        logger.warning(f'VCI FILE: {vci_file.filename}')
        logger.info(f'VCI FILE IS DIPLOID: {vci_file.is_diploid()}')
        vci_file.parse(reverse=False)

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

    total = 0
    success = 0
    fail = 0
    ret = OrderedDict()

    for r in all_regions:
        logger.debug(f'ORIGINAL: {str(r)}')

        total += 1

        if total % 10000 == 0:
            logger.info(f'Processed {total:,} regions')

        for lr in left_right:
            seq_id = f'{r.seq_id}{lr}'
            mappings = vci_file.find_mappings(
                seq_id, r.start - 1, r.end
            )

            # unmapped
            if mappings is None:
                ret[r] = None
                logger.debug(f'\t{r} -> No mappings')
                fail += 0
                continue
            else:
                logger.debug(f'\t{r} -> {len(mappings)} mappings')

            success += 1
            start = mappings[0].to_start + 1
            end = mappings[-1].to_end

            ret[r] = region.Region(seq_id, start, end)
            logger.debug(f'\t{r} -> {ret[r]}')

    logger.warning(f'Converted {success:,} of {total:,} records')

    return ret


#
# Collection of functions related to GTF
#

# standard library imports
from bz2 import BZ2File
from collections import namedtuple, OrderedDict
from dataclasses import dataclass
from gzip import GzipFile
from typing import TextIO
import sys
import time
import traceback
import pandas as pd

# 3rd party library imports
from intervaltree import IntervalTree

# local library imports
from g2gtools.exceptions import G2GGTFError
from g2gtools import g2g_utils
from g2gtools import vci


logger = g2g_utils.get_logger('g2gtools')


gtfInfoFields = [
    'seqid',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attributes',
]
GTFRecord = namedtuple('GTFRecord', gtfInfoFields)

ATTRIBUTES_TO_ALTER = {
    'gene_id': 'gene_id',
    'transcript_id': 'transcript_id',
    'exon_id': 'exon_id',
    'protein_id': 'protein_id',
    'ccds_id': 'ccds_id',
}

ATTRIBUTES_TO_ALTER_GFF = {
    'ID': 'ID',
    'Parent': 'Parent',
    'Name': 'Name',
    'Note': 'Note',
}


class GTF:
    """
    Simple GTF object for parsing GTF files.

    This class provides functionality to read and iterate through GTF
    (Gene Transfer Format) files, supporting transparent gzip decompression.
    It parses each line of a GTF file into a GTFRecord object with appropriate
    fields according to the GTF specification.

    Reference: http://blogger.nextgenetics.net/?e=27

    Attributes:
        file_name (str): Path to the GTF file.
        current_line (str): The current line being processed.
        current_record (GTFRecord | None): The current parsed GTF record.
        reader (GzipFile | TextIO | BZ2File | None): File reader object for the GTF file.
    """

    # Type annotations for instance variables
    file_name: str
    current_line: str
    current_record: GTFRecord | None
    reader: GzipFile | TextIO | BZ2File | None

    def __init__(self, file_name: str) -> None:
        """
        Initialize a new GTF file reader.

        Args:
            file_name: Path to the GTF file to read.
        """
        self.file_name = file_name
        self.current_line = ''
        self.current_record = None
        self.reader = g2g_utils.open_resource(file_name)

    def __iter__(self) -> 'GTF':
        """
        Make the GTF object iterable.

        Returns:
            The GTF object itself as an iterator.
        """
        return self

    def __next__(self) -> GTFRecord:
        """
        Get the next GTF record from the file.

        Reads the next line from the GTF file, skipping comment lines
        (those starting with '#' or '!'), parses it, and returns a
        GTFRecord object.

        Returns:
            A GTFRecord object representing the next record in the file.

        Raises:
            StopIteration: When the end of the file is reached.
        """
        self.current_line = g2g_utils.s(self.reader.__next__())

        # Skip comment lines
        while self.current_line.startswith(
            '#'
        ) or self.current_line.startswith('!'):
            self.current_line = g2g_utils.s(self.reader.__next__())

        # Parse the line into a GTFRecord
        self.current_record = parse_gtf_line(self.current_line)
        return self.current_record


def attributes_to_odict(attributes: str) -> OrderedDict[str, str]:
    """
    Parse the GTF attribute column and return a dictionary.

    Args:
        attributes: A string of attributes.

    Returns:
        A dictionary of attributes.
    """
    if attributes == '.':
        return OrderedDict()

    ret = OrderedDict()
    for attribute in attributes.strip().split('";'):
        if len(attribute):
            attribute = f'{attribute}"'
            elem = attribute.strip().split(' ')
            key = elem[0]
            val = ' '.join(elem[1:])

            if val[0] == '"':
                val = val[1:]
            if val[-1] == '"':
                val = val[0:-1]
            ret[key] = val
    return ret


def odict_to_attributes(attributes: dict[str, str]) -> str:
    """
    Parse the GTF attribute dictionary and return a string.

    Args:
        attributes: A dictionary of GTF attributes.

    Returns:
        A string representation.
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f'{k} "{v}"')
        temp_atts = '; '.join(atts)
        return temp_atts.rstrip() + ';'

    return '.'


def parse_gtf_line(line: str) -> GTFRecord:
    """
    Parse the GTF/GFF line.

    Args:
        line: A line from GTF file.

    Returns:
        A GTFRecord.
    """
    elem = line.strip().split('\t')

    # If this fails, the file format is not standard-compatible
    if len(elem) != len(gtfInfoFields):
        raise G2GGTFError('Improperly formatted GTF file')

    data = {
        'seqid': None if elem[0] == '.' else elem[0].replace('"', ''),
        'source': None if elem[1] == '.' else elem[1].replace('"', ''),
        'type': None if elem[2] == '.' else elem[2].replace('"', ''),
        'start': None if elem[3] == '.' else int(elem[3]),
        'end': None if elem[4] == '.' else int(elem[4]),
        'score': None if elem[5] == '.' else float(elem[5]),
        'strand': None if elem[6] == '.' else elem[6].replace('"', ''),
        'frame': None if elem[7] == '.' else elem[7].replace('"', ''),
    }

    # since GTF will contain '"' vs. GFF will contain '=' in the elem[8]
    if '"' in elem[8]:
        data.update({'attributes': attributes_to_odict(elem[8])})
    elif '=' in elem[8]:
        data.update({'attributes': attributes_to_odict_gff(elem[8])})

    return GTFRecord(**data)


def convert_gtf_file(
    vci_file: str | vci.VCIFile,
    gtf_file_name_in: str,
    gtf_file_name_out: str | None = None,
    reverse: bool | None = False,
) -> None:
    """
    Convert GTF file.

    Args:
        vci_file: Name of the VCI file or a VCIFile object.
        gtf_file_name_in: Input GTF file to convert.
        gtf_file_name_out: Name of output GTF file, None for stdout.
        reverse: True to process VCI in reverse.
    """
    start_time = time.time()

    if not isinstance(vci_file, vci.VCIFile):
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        vci_file.parse(reverse)

    logger.warning(f'Input VCI File: {vci_file.get_filename()}')
    gtf_file_name_in = g2g_utils.check_file(gtf_file_name_in)
    logger.warning(f'Input GTF File: {gtf_file_name_in}')
    logger.warning(f'Output GTF File: {gtf_file_name_out}')

    if gtf_file_name_out:
        output_file = g2g_utils.check_file(gtf_file_name_out, 'w')
        unmapped_file = f'{output_file}.unmapped'
        gtf_out = open(output_file, 'w')
        gtf_unmapped_file = open(unmapped_file, 'w')
        logger.info(f'Output GTF File: {output_file}')
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(gtf_file_name_in)
        unmapped_file = f'{input_name}.unmapped'
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        gtf_out = sys.stdout
        gtf_unmapped_file = open(unmapped_file, 'w')
        logger.info('Output GTF File: stdout')

    logger.info(f'Output UNMAPPED file: {unmapped_file}')
    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

    logger.info('Converting GTF file')

    gtf_file = GTF(gtf_file_name_in)

    total = 0
    success = 0
    fail = 0

    try:

        # GTF is 1 based, bx-python is 0 based
        # when we do the querying, we subtract 1 from the GTF file start position
        # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
        #       bx-python (s, e) does it s <= x < e.

        for record in gtf_file:
            logger.debug(f'\nORIGINAL: {str(gtf_file.current_line).strip()}')

            total += 1

            if total % 100000 == 0:
                logger.info(f'Processed {total:,} lines')

            for lr in left_right:
                seq_id = f'{record.seqid}{lr}'
                mappings = vci_file.find_mappings(
                    seq_id, record.start - 1, record.end
                )

                # unmapped
                if mappings is None:
                    logger.debug('Fail due to no mappings')
                    gtf_unmapped_file.write(gtf_file.current_line)
                    fail += 0
                    continue
                else:
                    logger.debug(f'{len(mappings)} mappings found')

                success += 1
                start = mappings[0].to_start + 1
                end = mappings[-1].to_end

                logger.debug(
                    f'({record.start - 1}, {record.end}) => ({start}, {end})'
                )

                elem = gtf_file.current_line.rstrip().split('\t')
                elem[0] = seq_id
                elem[3] = start
                elem[4] = end

                if lr:
                    attributes = attributes_to_odict(elem[8])
                    logger.debug(f'attributes={attributes}')
                    for k, v in ATTRIBUTES_TO_ALTER.items():
                        if k in attributes:
                            # only append '_L' or '_R' if it is not empty
                            if len(attributes[k]) > 0:
                                attributes[k] = f'{attributes[k]}{lr}'
                            else:
                                attributes[k] = f'{attributes[k]}'

                    elem[8] = odict_to_attributes(attributes)

                logger.debug('     NEW: {0}'.format('\t'.join(map(str, elem))))

                gtf_out.write('\t'.join(map(str, elem)))
                gtf_out.write('\n')

        if gtf_out:
            gtf_out.close()

        if gtf_unmapped_file:
            gtf_unmapped_file.close()

        if vci_file.is_haploid():
            logger.warning(
                f'Converted {success:,} records from {total:,} records'
            )
        else:
            logger.warning(
                f'Converted {int(success/2):,} records from {total:,} records'
            )

        logger.warning('GTF file converted')
    except KeyboardInterrupt:
        logger.warning('Keyboard interrupt detected. Exiting...')
    finally:
        fmt_time = g2g_utils.format_time(start_time, time.time())
        logger.warning(f'Time: {fmt_time}')


def convert_gff_file(
    vci_file: str | vci.VCIFile,
    gff_file_name_in: str,
    gff_file_name_out: str | None = None,
    reverse: bool | None = False,
) -> None:
    """
    Convert GFF file.

    Args:
        vci_file: Name of the VCI file or a VCIFile object.
        gff_file_name_in: Input GFF file to convert.
        gff_file_name_out: Name of output GFF file, None for stdout.
        reverse: True to process VCI in reverse.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start_time = time.time()

    if not isinstance(vci_file, vci.VCIFile):
        vci_file = g2g_utils.check_file(vci_file)
        vci_file = vci.VCIFile(vci_file)
        vci_file.parse(reverse)

    logger.warning(f'Input VCI File: {vci_file.get_filename()}')
    gff_file_name_in = g2g_utils.check_file(gff_file_name_in)
    logger.warning(f'Input GTF File: {gff_file_name_in}')
    logger.warning(f'Output GTF File: {gff_file_name_out}')

    if gff_file_name_out:
        output_file = g2g_utils.check_file(gff_file_name_out, 'w')
        unmapped_file = f'{output_file}.unmapped'
        gff_out = open(output_file, 'w')
        gff_unmapped_file = open(unmapped_file, 'w')
        logger.info(f'Output GTF File: {output_file}')
    else:
        input_dir, input_name = g2g_utils.get_dir_and_file(gff_file_name_in)
        unmapped_file = f'{input_name}.unmapped'
        unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
        gff_out = sys.stdout
        gff_unmapped_file = open(unmapped_file, 'w')
        logger.info('Output GTF File: stdout')

    logger.info(f'Output UNMAPPED file: {unmapped_file}')

    left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

    logger.info('Converting GFF file')

    gff_file = GTF(gff_file_name_in)

    total = 0
    success = 0
    fail = 0

    try:
        # GTF is 1 based, bx-python is 0 based
        # when we do the querying, we subtract 1 from the GTF file start position
        # K.B.  Also note in gtf when (s, e) is given...it should mean s <= x <= e.
        #       bx-python (s, e) does it s <= x < e.
        for record in gff_file:
            logger.debug(f'ORIGINAL: {str(gff_file.current_line).strip()}')

            total += 1

            if total % 100000 == 0:
                logger.info(f'Processed {total:,} lines')

            for lr in left_right:
                seq_id = f'{record.seqid}{lr}'
                mappings = vci_file.find_mappings(
                    seq_id, record.start - 1, record.end
                )

                # unmapped
                if mappings is None:
                    logger.debug('\tFail due to no mappings')
                    gff_unmapped_file.write(gff_file.current_line)
                    fail += 0
                    continue
                else:
                    logger.debug(f'{len(mappings)} mappings found')

                success += 1
                start = mappings[0].to_start + 1
                end = mappings[-1].to_end

                logger.debug(
                    f'({record.start}, {record.end}) => ({start}, {end})'
                )

                elem = gff_file.current_line.rstrip().split('\t')
                elem[0] = seq_id
                elem[3] = start
                elem[4] = end

                if lr:
                    attributes = attributes_to_odict_gff(elem[8])
                    for k, v in ATTRIBUTES_TO_ALTER_GFF.items():
                        if k in attributes:
                            attributes[k] = f'{attributes[k]}{lr}'
                    elem[8] = odict_to_attributes_gff(attributes)

                logger.debug('     NEW: {0}'.format('\t'.join(map(str, elem))))

                gff_out.write('\t'.join(map(str, elem)))
                gff_out.write('\n')

                if gff_out:
                    gff_out.close()

                if gff_unmapped_file:
                    gff_unmapped_file.close()

                if vci_file.is_haploid():
                    logger.warning(
                        f'Converted {success:,} records from {total:,} records'
                    )
                else:
                    logger.warning(
                        f'Converted {int(success / 2):,} records from {total:,} records'
                    )

                logger.warning('GFF file converted')
    except KeyboardInterrupt:
        logger.warning('Keyboard interrupt detected. Exiting...')
    finally:
        fmt_time = g2g_utils.format_time(start_time, time.time())
        logger.warning(f'Time: {fmt_time}')


def attributes_to_odict_gff(attributes: str) -> OrderedDict[str, str]:
    """
    Parse the GFF attribute column and return a dict.

    Args:
        attributes: The string of attributes.

    Returns:
        A dictionary of attributes.
    """
    if attributes == '.':
        return OrderedDict()

    ret_gff = OrderedDict()
    for attribute in attributes.strip().split(';'):
        if len(attribute):
            elem = attribute.strip().split('=')
            key = elem[0]
            val = ','.join(elem[1:])

            ret_gff[key] = val

    return ret_gff


def odict_to_attributes_gff(attributes: dict[str, str]) -> str:
    """
    Assemble the GFF attributes into a string.

    Args:
        attributes: A dictionary of attributes.

    Returns:
        A string version of the GFF attributes.
    """
    if attributes:
        atts = []
        for k, v in attributes.items():
            atts.append(f'{k}={v}')
        temp_atts = ';'.join(atts)
        return temp_atts.rstrip()

    return '.'


@dataclass
class GeneInfo:
    """
    Store gene information from GTF file.

    This class represents a gene entry parsed from a GTF (Gene Transfer Format) file,
    containing genomic coordinates and identifiers for a specific gene. It also
    tracks variant counts within the gene region.

    Attributes:
        gene_id (str): The unique identifier for the gene (e.g., "ENSG00000139618").
        gene_name (str): The common name or symbol of the gene (e.g., "BRCA2").
        chrom (str): The chromosome name where this gene is located (e.g., "chr13").
        start (int): The starting position of the gene on the chromosome.
        end (int): The ending position of the gene on the chromosome (inclusive).
        snp_count (int): Number of single nucleotide polymorphisms (SNPs) in this gene.
                         Defaults to 0.
        indel_count (int): Number of insertions and deletions (indels) in this gene.
                           Defaults to 0.
    """

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    snp_count: int = 0
    indel_count: int = 0


@dataclass
class TranscriptInfo:
    """
    Store transcript information from GTF file.

    This class represents a transcript entry parsed from a GTF (Gene Transfer Format) file,
    containing genomic coordinates and identifiers for a specific transcript. It also
    tracks variant counts within the transcript region.

    Attributes:
        gene_id (str): The unique identifier for the gene this transcript belongs to.
        gene_name (str): The common name or symbol of the gene (e.g., "BRCA1").
        transcript_id (str): The unique identifier for this specific transcript.
        chrom (str): The chromosome name where this transcript is located (e.g., "chr1").
        start (int): The starting position of the transcript on the chromosome.
        end (int): The ending position of the transcript on the chromosome (inclusive).
        snp_count (int): Number of single nucleotide polymorphisms (SNPs) in this transcript.
                         Defaults to 0.
        indel_count (int): Number of insertions and deletions (indels) in this transcript.
                           Defaults to 0.
    """

    gene_id: str
    gene_name: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    snp_count: int = 0
    indel_count: int = 0


def parse_gtf(gtf_file: str, chrom: str = None):
    """
    Parse GTF file to extract gene information.

    Args:
        gtf_file: Path to the GTF file

    Returns:
        Dictionary mapping Ensembl gene IDs to GeneInfo objects
    """
    logger.info(f'Parsing GTF file: {gtf_file}')
    transcripts = {}
    transcripts_tree = IntervalTree()

    try:
        # Read GTF file with pandas for faster processing
        # Only read the columns we need
        cols = [0, 2, 3, 4, 8]  # seqname, feature, start, end, attributes
        df = pd.read_csv(
            gtf_file,
            sep='\t',
            comment='#',
            header=None,
            usecols=cols,
            names=['chrom', 'feature', 'start', 'end', 'attributes'],
            low_memory=False,  # Avoid DtypeWarning
        )
        df['chrom'] = df['chrom'].astype(
            str
        )  # Ensure chromosome is string type

        # Filter for gene entries only
        # filtered_df = df[df['feature'] == 'gene']
        filtered_df = df[df['feature'] == 'transcript']

        logger.debug(f'Found {len(filtered_df)} entries in GTF file')

        # Process each gene entry
        for _, row in filtered_df.iterrows():
            attributes = row['attributes']

            # Extract gene_id and gene_name from attributes
            gene_id_match = None
            gene_name_match = None
            transcript_id_match = None

            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id_match = attr.split('"')[1]
                elif attr.startswith('transcript_id'):
                    transcript_id_match = attr.split('"')[1]
                elif attr.startswith('gene_name'):
                    gene_name_match = attr.split('"')[1]

            if gene_id_match and gene_name_match and transcript_id_match:
                if (chrom is None) or (chrom and row['chrom'] == chrom):

                    # Add gene to dictionary
                    # GTF is 1-based, inclusive
                    transcripts[transcript_id_match] = TranscriptInfo(
                        gene_id=gene_id_match,
                        gene_name=gene_name_match,
                        transcript_id=transcript_id_match,
                        chrom=row['chrom'],
                        start=int(row['start']),
                        end=int(row['end']),
                        snp_count=0,
                        indel_count=0,
                    )

                    transcripts_tree[
                        int(row['start']) : int(row['end'])
                    ] = transcript_id_match

        logger.info(
            f'Successfully parsed {len(transcripts)} transcripts from GTF file'
        )

        return {'features': transcripts, 'tree': transcripts_tree}

    except Exception as e:
        logger.error(f'Error parsing GTF file: {e}')
        logger.error(traceback.format_exc())
        sys.exit(1)

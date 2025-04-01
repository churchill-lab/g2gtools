"""
Collection of functions related to FASTA files.
"""

# standard library imports
from io import StringIO
from pathlib import Path
from typing import Any, Iterator
import collections
import difflib
import os
import re
import sys
import time

# 3rd party library imports
import pysam

# local library imports
from g2gtools.exceptions import G2GFastaError
from g2gtools.exceptions import G2GRegionError
from g2gtools.exceptions import G2GValueError
import g2gtools.gtf_db as gtf_db
import g2gtools.region as region
import g2gtools.g2g_utils as g2g_utils


logger = g2g_utils.get_logger('g2gtools')


REGEX_FASTA_HEADER = re.compile(r'>(\S*)\s*(.*)')
FASTA_HEADER_FIELDS = ['id', 'description']
FastaHeader = collections.namedtuple('FastaHeader', FASTA_HEADER_FIELDS)


class FAIEntry:
    """
    Represents an entry in a FASTA index file.

    This class stores information about a single sequence in a FASTA file,
    including its identifier, length, and file offset information needed
    for random access.

    Attributes:
        idx (str | None): Sequence identifier.
        length (int | None): Length of the sequence.
        offset (int | None): Byte offset of the sequence in the FASTA file.
        line_len (int | None): Number of bases per line.
        line_len_bytes (int | None): Number of bytes per line including newlines.
    """
    idx: str | None
    length: int | None
    offset: int | None
    line_len: int | None
    line_len_bytes: int | None

    def __init__(
        self,
        idx: str | None = None,
        length: int | None = 0,
        offset: int | None = -1,
        line_len: int | None = None,
        line_len_bytes: int | None = None,
    ) -> None:
        """
        Initialize a FAIEntry object.

        Args:
            idx: The sequence identifier.
            length: The length of the sequence.
            offset: The byte offset of the sequence in the FASTA file.
            line_len: The number of bases per line.
            line_len_bytes: The number of bytes per line including newlines.
        """
        self.idx = idx
        self.length = length
        self.offset = offset
        self.line_len = line_len
        self.line_len_bytes = line_len_bytes

    def __str__(self) -> str:
        """
        Get a string representation of the FAIEntry.

        Returns:
            A tab-delimited string representation of the entry.
        """
        return (
            f'{self.idx}\t{self.length}\t{self.offset}\t'
            f'{self.line_len}\t{self.line_len_bytes}'
        )


class FAI:
    """
    Represents a FASTA index file (.fai).

    This class provides access to the index of a FASTA file, which enables
    efficient random access to sequences within the file. It reads and parses
    the .fai file format used by samtools and pysam.

    Attributes:
        fasta_file (str): Path to the FASTA file.
        fai_file (str): Path to the FASTA index file.
        records (collections.OrderedDict): Dictionary mapping sequence IDs to
        FAIEntry objects.
    """
    fasta_file: str
    fai_file: str
    records: collections.OrderedDict[str, FAIEntry]

    def __init__(self, fasta_file_name: str) -> None:
        """
        Initialize a FASTA index object.

        Reads the .fai file associated with the given FASTA file. Tries both
        .fai and .gz.fai extensions.

        Args:
            fasta_file_name: The path to the FASTA file.

        Raises:
            G2GFastaError: If the index file cannot be found.
        """
        self.fasta_file = fasta_file_name
        self.fai_file = f'{fasta_file_name}.fai'
        self.records = collections.OrderedDict()

        if os.path.exists(self.fai_file):
            self.read()
        else:
            self.fai_file = f'{fasta_file_name}.gz.fai'
            if os.path.exists(self.fai_file):
                self.read()
            else:
                raise G2GFastaError('Unable to find fasta index file')

    def __getitem__(self, entry: str | int) -> FAIEntry:
        """
        Get a FAIEntry by sequence ID or index.

        Args:
            entry: Sequence ID (str) or index position (int).

        Returns:
            The FAIEntry object for the specified sequence.

        Raises:
            G2GFastaError: If the entry is not found in the index.
        """
        try:
            if isinstance(entry, int):
                entry = tuple(self.records.keys())[entry]
            return self.records[entry]
        except KeyError:
            raise G2GFastaError(f'{entry} not in {self.fai_file}.')

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over sequence IDs in the index.

        Returns:
            An iterator over sequence IDs.
        """
        for keys in self.records:
            yield keys

    def __repr__(self) -> str:
        """
        Get a string representation of the FAI object.

        Returns:
            A string representation including the class name and FAI file path.
        """
        return f'{self.__class__}("{self.fai_file}")'

    def read(self) -> None:
        """
        Read and parse the FASTA index file.

        Populates the records dictionary with FAIEntry objects.

        Raises:
            ValueError: If a duplicate sequence ID is found in the index.
        """
        with open(self.fai_file) as index:
            for line in index:
                line = line.strip()
                idx, length, offset, line_len, line_len_bytes = line.split(
                    '\t'
                )
                if idx in self.records:
                    raise ValueError(f'A duplicate entry was found: {idx}')
                else:
                    self.records[idx] = FAIEntry(
                        idx,
                        int(length),
                        int(offset),
                        int(line_len),
                        int(line_len_bytes),
                    )

    def get_pos(
        self, seq_id: str, start: int, end: int
    ) -> tuple[int, int, int]:
        """
        Get the byte positions in the FASTA file for a specified sequence region.

        Calculates the byte offsets needed for random access to a specific
        region of a sequence, accounting for line breaks in the FASTA format.

        Args:
            seq_id: The sequence identifier.
            start: The start position (0-based).
            end: The end position (exclusive).

        Returns:
            A tuple of (starting_byte, ending_byte, byte_length).
        """
        chrom = self.records[seq_id]
        fai_entry_length = chrom.length
        fai_entry_offset = chrom.offset
        fai_entry_line_length = chrom.line_len
        fai_entry_line_length_bytes = chrom.line_len_bytes

        seq_len = end - start
        fai_diff = fai_entry_line_length_bytes - fai_entry_line_length
        line_ratio = fai_entry_line_length * fai_diff
        newlines_total = int(fai_entry_length / line_ratio)

        newlines_before = 0
        if start > 0:
            newlines_before = int(start / line_ratio)

        newlines_to_end = int(end / line_ratio)
        byte_len_seq = newlines_to_end - newlines_before + seq_len
        byte_start = fai_entry_offset + newlines_before + start
        byte_end = fai_entry_offset + newlines_total + fai_entry_length

        return byte_start, byte_end, byte_len_seq

    def dump(self) -> None:
        """
        Print a human-readable representation of the FASTA index.

        Displays information about the FASTA file, index file, and all records.
        """
        print('FAI')
        print(f'Fasta File: {self.fasta_file}')
        print(f'FAI File: {self.fai_file}')
        for record in self.records:
            print(str(self.records[record]))


class FastaFile:
    """
    Encapsulate a Fasta file.

    This class provides a wrapper around pysam.FastaFile with additional
    functionality for handling haploid and diploid genomes. It uses delegation
    to provide access to the underlying pysam.FastaFile methods.

    Attributes:
        filename (str): Path to the FASTA file.
        _fasta_file (pysam.FastaFile): The underlying pysam FASTA file object.
        fai (FAI): The FASTA index object.
        haploid (bool): Whether the FASTA file represents a haploid genome.
        diploid (bool): Whether the FASTA file represents a diploid genome.
    """
    filename: str
    _fasta_file: pysam.FastaFile
    fai: FAI
    haploid: bool
    diploid: bool

    def __init__(self, filename: str) -> None:
        """
        Initialize a FastaFile object.

        Args:
            filename: Path to the FASTA file.
        """
        self.filename = filename
        self._fasta_file = pysam.FastaFile(self.filename)
        self.fai = FAI(self.filename)
        self.haploid = False
        self.diploid = False

        # determine if the FASTA file is haploid or diploid
        for record in self.fai.records:
            if record.endswith(('_L', '_R')):
                self.diploid = True
                self.haploid = False
                break
        if not self.diploid:
            self.haploid = True

    def get_filename(self) -> str:
        """
        Get the filename of the Fasta File.

        Returns:
            The filename.
        """
        return self.filename

    def is_diploid(self) -> bool:
        """
        Check if this FASTA file represents a diploid genome.

        Returns:
            True if this is a diploid genome, False otherwise.
        """
        return self.diploid

    def is_haploid(self) -> bool:
        """
        Check if this FASTA file represents a haploid genome.

        Returns:
            True if this is a haploid genome, False otherwise.
        """
        return self.haploid

    def fetch_list(
        self,
        reference: str | None = None,
        start: int | None = None,
        end: int | None = None,
        region: str | None = None,
    ) -> list[str]:
        """
        Fetch sequence as a list of characters.

        Args:
            reference: The reference sequence name.
            start: The start position (0-based).
            end: The end position (exclusive).
            region: Region string in format "chr:start-end".

        Returns:
            A list of characters representing the sequence.
        """
        _start = 0 if start < 0 else start
        return list(
            self._fasta_file.fetch(
                reference=reference, start=_start, end=end, region=region
            )
        )

    def __getattr__(self, name: str) -> Any:
        """
        Delegate attribute access to the underlying pysam.FastaFile object.

        Args:
            name: The attribute name to access.

        Returns:
            The attribute value from the underlying pysam.FastaFile object.
        """
        return getattr(self._fasta_file, name)


def extract(
    fasta_file: str | FastaFile,
    locations: region.Region | list[region.Region],
    output_file_name: str | None = None,
    reverse: bool | None = False,
    complement: bool | None = False,
    raw: bool | None = False,
):
    """
    Extract the Fasta sequences.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        locations: Either a list of Regions or just a Region.
        output_file_name: Name of output file, None for sys.stdout.
        reverse: True to reverse the extracted sequence.
        complement: True to complement the extracted sequence.
        raw: True for just the sequence to be output, no fasta line.

    Raises:
        G2GRegionError: When a Region is not correct.
    """
    start_time = time.time()

    if not isinstance(fasta_file, pysam.FastaFile):
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta_file = pysam.FastaFile(fasta_file)

    logger.warning(
        f'Input Fasta File: {g2g_utils.convert_if_bytes(fasta_file.filename)}'
    )

    if output_file_name:
        output_file_name = g2g_utils.check_file(output_file_name, 'w')
        logger.warning(f'Output Fasta File: {output_file_name}')
        fasta_out = open(output_file_name, 'w')
    else:
        logger.warning('Output Fasta File: stderr')
        fasta_out = sys.stderr

    num_sequences = 0

    try:
        if not isinstance(locations, list):
            locations = [locations]

        for location in locations:
            logger.debug(f'LOCATION: {location}')
            if not isinstance(location, region.Region):
                raise G2GRegionError(
                    f'Found: {location}, which is not a Region'
                )

            if location.seq_id not in fasta_file.references:
                continue

            sequence = fasta_file.fetch(
                location.seq_id, location.start, location.end
            )

            if location.strand == '+':
                if reverse and not complement:
                    sequence = g2g_utils.reverse_sequence(sequence)
                elif not reverse and complement:
                    sequence = g2g_utils.complement_sequence(sequence)
                elif reverse and complement:
                    sequence = g2g_utils.reverse_complement_sequence(sequence)
            else:
                sequence = g2g_utils.reverse_complement_sequence(sequence)

            if raw:
                fasta_out.write(sequence)
                fasta_out.write('\n')
            else:
                # internally storing 0 base, output 1 base
                start = location.start + 1
                end = location.end
                len_seq = len(sequence)

                # if someone requests something to long, just truncate
                # the length
                if len_seq < end:
                    end = len_seq + start - 1

                fasta_id = f'>{location.seq_id}:{start}-{end}'
                logger.debug(f'Fasta ID: {fasta_id}')
                logger.debug(f'{sequence[:10]}...{sequence[-10:]}')
                num_sequences += 1

                if location.name:
                    fasta_id = (
                        f'>{location.name} {location.seq_id}:{start}-{end}'
                    )
                    logger.debug(f'Location name used, Fasta ID: {fasta_id}')

                fasta_out.write(f'{fasta_id}\n')
                g2g_utils.write_sequence(sequence, fasta_out)

        logger.warning(f'Extracted {num_sequences:,} sequences')
        logger.warning(f'Fasta file created')

    except G2GRegionError as e:
        logger.info(str(e).rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(str(e).rstrip())
        raise e
    finally:
        fasta_out.close()
        fmt_time = g2g_utils.format_time(start_time, time.time())
        logger.warning(f'Time: {fmt_time}')


def extract_id(
    fasta_file: str | FastaFile,
    identifier: str,
    output_file_name: str | None = None,
    reverse: bool | None = False,
    complement: bool | None = False,
    raw: bool | None = False,
):
    """
    Extract the Fasta sequences for a given identifier.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        identifier: An identifier.
        output_file_name: Name of output file, None for sys.stdout.
        reverse: True to reverse the extracted sequence.
        complement: True to complement the extracted sequence.
        raw: True for just the sequence to be output, no fasta line.
    """
    start_time = time.time()

    if not isinstance(fasta_file, pysam.FastaFile):
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta_file = pysam.FastaFile(fasta_file)

    logger.warning(
        f'Input Fasta File: {g2g_utils.convert_if_bytes(fasta_file.filename)}'
    )

    if output_file_name:
        output_file_name = g2g_utils.check_file(output_file_name, 'w')
        logger.warning(f'Output Fasta File: {output_file_name}')
        fasta_out = open(output_file_name, 'w')
    else:
        logger.warning('Output Fasta File: stderr')
        fasta_out = sys.stderr

    try:
        sequence = fasta_file[identifier]

        if reverse and not complement:
            sequence = g2g_utils.reverse_sequence(sequence)
        elif not reverse and complement:
            sequence = g2g_utils.complement_sequence(sequence)
        elif reverse and complement:
            sequence = g2g_utils.reverse_complement_sequence(sequence)

        if raw:
            fasta_out.write(sequence)
        else:
            fasta_id = f'>{identifier}\n'
            fasta_out.write(fasta_id)

            for line in g2g_utils.wrap_sequence(sequence):
                fasta_out.write(line.strip())
                fasta_out.write('\n')

    except G2GRegionError as e:
        logger.info(str(e).rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(str(e).rstrip())
        raise e
    finally:
        if output_file_name:
            fasta_out.close()
        fmt_time = g2g_utils.format_time(start_time, time.time())
        logger.warning(f'Time: {fmt_time}')




def diff_files(file_name_1: str, file_name_2: str, verbose: bool = True) -> tuple[set[str], set[str], set[str]]:
    """
    Compare two FASTA files and identify sequences that are identical, different, or unique to each file.

    This function opens two FASTA files and compares their sequences by reference ID.
    It identifies which sequences are identical between the files, which have the same ID
    but different content, and which are unique to each file.

    Args:
        file_name_1: Path to the first FASTA file
        file_name_2: Path to the second FASTA file
        verbose: Whether to print summary statistics to stdout (default: True)

    Returns:
        A tuple containing three sets:
        - Set of sequence IDs that are identical between the files
        - Set of sequence IDs that have different content between the files
        - Set of sequence IDs that are unique to either file

    Raises:
        FileNotFoundError: If either of the specified files doesn't exist
        ValueError: If the files cannot be opened as FASTA files
    """
    # Validate input files
    for file_path in (file_name_1, file_name_2):
        if not Path(file_path).exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

    # Open FASTA files using context managers for proper resource handling
    with pysam.FastaFile(file_name_1) as fasta_1, pysam.FastaFile(file_name_2) as fasta_2:
        # Get all reference IDs from both files
        refs_1 = set(fasta_1.references)
        refs_2 = set(fasta_2.references)

        # Find common references and unique references
        common_refs = refs_1.intersection(refs_2)
        unique_refs = refs_1.symmetric_difference(refs_2)

        # Compare sequences for common references
        identical_seqs = set()
        different_seqs = set()

        for seq_id in common_refs:
            if fasta_1.fetch(seq_id) == fasta_2.fetch(seq_id):
                identical_seqs.add(seq_id)
            else:
                different_seqs.add(seq_id)

    # Print summary if verbose mode is enabled
    if verbose:
        print(f"# EQUAL: {len(identical_seqs)}")
        print(f"# DIFF: {len(different_seqs)}")
        print(f"# UNIQUE: {len(unique_refs)}")

        if different_seqs:
            print("\nSequences with different content:")
            for seq in sorted(different_seqs):
                print(f"  {seq}")

        if unique_refs:
            print("\nUnique sequences:")
            for seq in sorted(unique_refs):
                source = "file 1" if seq in refs_1 else "file 2"
                print(f"  {seq} (only in {source})")

    return identical_seqs, different_seqs, unique_refs


def diff_sequence(sequence_1: str, sequence_2: str) -> None:
    """
    Compare two sequences.

    Example:
        cases=[('afrykanerskojęzyczny', 'afrykanerskojęzycznym'),
               ('afrykanerskojęzyczni', 'nieafrykanerskojęzyczni'),
               ('afrykanerskojęzycznym', 'afrykanerskojęzyczny'),
               ('nieafrykanerskojęzyczni', 'afrykanerskojęzyczni'),
               ('nieafrynerskojęzyczni', 'afrykanerskojzyczni'),
               ('abcdefg','xac')]

        for a,b in cases:
            print('{} => {}'.format(a,b))
            for i,s in enumerate(difflib.ndiff(a, b)):
                if s[0]==' ': continue
                elif s[0]=='-':
                    print(u'Delete "{}" from position {}'.format(s[-1],i))
                elif s[0]=='+':
                    print(u'Add "{}" to position {}'.format(s[-1],i))
            print()

    Args:
        sequence_1: The first sequence.
        sequence_2: The second sequence.
    """
    for i, s in enumerate(difflib.ndiff(sequence_1, sequence_2)):
        print(i, s)

        if s[0] == ' ':
            continue
        elif s[0] == '-':
            print(f'Delete "{s[-1]}" from position {i}')
        elif s[0] == '+':
            print(f'Add "{s[-1]}" to position {i}')
    print()


def get_pos(
    fai: FAI, chrom: str, start: int, end: int
) -> tuple[int, int, int]:
    """
    Get the byte positions of the fasta file using the specified location

    Args:
        fai: The FAI object (Fasta file index).
        chrom: The chromosome.
        start: The start position.
        end: The end position.

    Returns:
        A tuple of starting byte, ending byte, and byte length.
    """
    chrom = fai.records[chrom]
    fai_entry_length = chrom.length
    fai_entry_offset = chrom.offset
    fai_entry_line_length = chrom.line_len
    fai_entry_line_length_bytes = chrom.line_len_bytes
    seq_len = end - start
    fai_diff = fai_entry_line_length_bytes - fai_entry_line_length
    line_ratio = fai_entry_line_length * fai_diff
    newlines_total = int(fai_entry_length / line_ratio)
    newlines_before = 0

    if start > 0:
        newlines_before = int(start / line_ratio)

    newlines_to_end = int(end / line_ratio)
    byte_len_seq = newlines_to_end - newlines_before + seq_len
    byte_start = fai_entry_offset + newlines_before + start
    byte_end = fai_entry_offset + newlines_total + fai_entry_length

    return byte_start, byte_end, byte_len_seq


def reformat(
    fasta_file_name: str, output_file_name: str = None, length: int = 60
) -> None:
    """
    Reformat a Fasta file to specified maximum line length.

    Args:
        fasta_file_name: The name of the Fasta file.
        output_file_name: Name of the output file or None for stdout.
        length: Maximum line length.
    """
    fasta_file_name = g2g_utils.check_file(fasta_file_name)
    output = sys.stdout

    if output_file_name:
        output_file = g2g_utils.check_file(output_file_name, 'w')
        output = open(output_file, 'w')

    with open(fasta_file_name, 'r') as fasta_fd:
        new_sequence = StringIO()
        new_header = '>UNKNOWN'

        for line in fasta_fd:
            line = line.strip()

            if line[0] == '>':
                if len(new_sequence.getvalue()) > 0:
                    output.write(new_header)
                    output.write('\n')
                    g2g_utils.write_sequence(
                        new_sequence.getvalue(), output, length
                    )

                    logger.info(f'Reformatting {new_header}')

                new_sequence = StringIO()
                new_header = line
            else:
                new_sequence.write(line)

        logger.info(f'Reformatting {new_header}')
        output.write(new_header)
        output.write('\n')
        g2g_utils.write_sequence(new_sequence.getvalue(), output, length)

    output.close()


def fasta_extract_transcripts(
    fasta_file: str | FastaFile,
    database_file_name: str,
    output_file_name: str | None = None,
    raw: bool | None = False,
) -> None:
    """
    Extract the transcripts sequences from the fasta_file given the database.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        database_file_name: Name of the database file.
        output_file_name: Name of the output file or None for stdout.
        raw: True to omit Fasta sequence name and/or info.
    """
    start = time.time()

    if not isinstance(fasta_file, FastaFile):
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta_file = FastaFile(fasta_file)

    logger.warning(f'Input Fasta File: {fasta_file.get_filename()}')

    database_file_name = g2g_utils.check_file(database_file_name)

    if output_file_name:
        output_file_name = g2g_utils.check_file(output_file_name, 'w')
        logger.warning(f'Output Fasta File: {output_file_name}')
        fasta_out = open(output_file_name, 'w')
    else:
        logger.warning('Output Fasta File: stderr')
        fasta_out = sys.stderr

    num_sequences = 0

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file_name)

        for i, transcript in enumerate(transcripts):
            logger.debug(f'Transcript={transcript}')

            if transcript.seqid not in fasta_file.references:
                logger.debug(
                    'Skipping, transcript seqid not in fasta references'
                )
                logger.debug(
                    f'{transcript.seqid} not in {fasta_file.references}'
                )
                continue

            new_sequence = StringIO()
            locations = []
            num_sequences += 1
            for ensembl_id, exon in transcript.exons.items():
                logger.debug(f'Exon ID={ensembl_id};{exon}')

                partial_seq = fasta_file.fetch(
                    exon.seqid, exon.start - 1, exon.end
                )
                partial_seq_str = str(partial_seq)

                if transcript.strand == 1:
                    new_sequence.write(partial_seq_str)
                else:
                    partial_seq_str = str(
                        g2g_utils.reverse_complement_sequence(partial_seq)
                    )
                    new_sequence.write(partial_seq_str)

                logger.debug(
                    f'{exon.seqid}:{exon.start}-{exon.end} '
                    f'(Length: {len(partial_seq)})\n{partial_seq_str}'
                )
                locations.append(f'{exon.seqid}:{exon.start}-{exon.end}')

            if raw:
                fasta_out.write(new_sequence.getvalue())
            else:
                t_strand = '-' if transcript.strand == -1 else '+'
                loc = '|'.join(locations)
                fasta_id = f'>{transcript.ensembl_id} {t_strand}|{loc}\n'
                fasta_out.write(fasta_id)

                for line in g2g_utils.wrap_sequence(new_sequence.getvalue()):
                    fasta_out.write(line.strip())
                    fasta_out.write('\n')

        logger.warning(f'Extracted {num_sequences:,} sequences')
        logger.warning(f'Fasta file created')

    except G2GValueError as e:
        logger.info(str(e).rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(str(e).rstrip())
        raise e
    finally:
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')


def fasta_extract_exons(
    fasta_file: str | FastaFile,
    database_file_name: str,
    output: str | None,
    raw: bool | None = False,
) -> None:
    """
    Extract the exons sequences from the fasta_file given the database.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        database_file_name: Name of the database file.
        output: Name of the output file or None for stdout.
        raw: True to omit Fasta sequence name and/or info.
    """
    start = time.time()

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file_name = g2g_utils.check_file(database_file_name)
    fasta_out = sys.stdout

    if output:
        output = g2g_utils.check_file(output, 'w')
        fasta_out = open(output, 'w')

    logger.warning(f'Input Fasta File: {fasta.filename}')
    logger.warning(f'Input DB File: {database_file_name}')
    logger.warning(f'Output Fasta File: {fasta_out.name}')
    num_exons = 0

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file_name)
        for i, transcript in enumerate(transcripts):
            if transcript.seqid not in fasta.references:
                continue

            for ensembl_id, exon in transcript.exons.items():
                logger.debug(f'Exon={exon}')
                num_exons += 1

                partial_seq = fasta.fetch(exon.seqid, exon.start - 1, exon.end)
                partial_seq_str = partial_seq

                if transcript.strand == -1:
                    rev_seq = g2g_utils.reverse_complement_sequence(
                        partial_seq
                    )
                    partial_seq_str = str(rev_seq)

                logger.debug(
                    f'{exon.seqid}:{exon.start}-{exon.end} '
                    f'(Length: {len(partial_seq)})\n{partial_seq_str}'
                )

                if raw:
                    fasta_out.write(partial_seq_str)
                else:
                    fasta_id = (
                        f'>{exon.ensembl_id} '
                        f'{exon.seqid}:{exon.start}-{exon.end}\n'
                    )
                    fasta_out.write(fasta_id)

                    for line in g2g_utils.wrap_sequence(partial_seq_str):
                        fasta_out.write(line.strip())
                        fasta_out.write('\n')

        logger.warning(f'Extracted {num_exons:,} exons')
    except G2GValueError as e:
        logger.info(str(e).rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(str(e).rstrip())
        raise e
    finally:
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')

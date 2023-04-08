"""
Collection of functions related to FASTA files.
"""
# standard library imports
from io import StringIO
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
import g2gtools.g2g as g2g
import g2gtools.g2g_utils as g2g_utils


REGEX_FASTA_HEADER = re.compile(r">(\S*)\s*(.*)")
FASTA_HEADER_FIELDS = ["id", "description"]
FastaHeader = collections.namedtuple("FastaHeader", FASTA_HEADER_FIELDS)


class FastaFile(object):
    """
    Encapsulate a Fasta file

    Delegation used with `pysam.FastaFile`
    """

    def __init__(self, filename):
        self.filename = filename
        self._fasta_file = None
        self._fasta_file = pysam.FastaFile(self.filename)

        self.fai = FAI(self.filename)
        self.haploid = False
        self.diploid = False

        for record in self.fai.records:
            if record.endswith(("_L", "_R")):
                self.diploid = True
                self.haploid = False
                break

        if not self.diploid:
            self.haploid = True

    def is_diploid(self):
        return self.diploid

    def is_haploid(self):
        return self.haploid

    def fetch_list(self, reference=None, start=None, end=None, region=None):
        _start = 0 if start < 0 else start
        return list(
            self._fasta_file.fetch(
                reference=reference, start=_start, end=end, region=region
            )
        )

    def __getattr__(self, name):
        return getattr(self._fasta_file, name)


class FAIEntry(object):
    def __init__(
        self,
        idx: str | None = None,
        length: int | None = 0,
        offset: int | None = -1,
        line_len: int | None = None,
        line_len_bytes: int | None = None,
    ):
        """
        Initialize a FAIEntry object.

        Args:
            idx: The index.
            length: The length.
            offset: The entry offset.
            line_len: The length of the line.
            line_len_bytes: The length of the line in bytes.
        """
        self.idx = idx
        self.length = length
        self.offset = offset
        self.line_len = line_len
        self.line_len_bytes = line_len_bytes

    def __str__(self):
        return (
            f"{self.idx}\t{self.length}\t{self.offset}\t"
            f"{self.line_len}\t{self.line_len_bytes}"
        )


class FAI(object):
    def __init__(self, fasta_file_name: str):
        """
        Initialize a `.fasta.FAI` file, an index into the Fasta file.

        Args:
            fasta_file_name: The name of the Fasta file.

        Returns:
            The initialized FAI object.
        """
        self.fasta_file = fasta_file_name
        self.fai_file = f"{fasta_file_name}.fai"
        self.records = collections.OrderedDict()

        if os.path.exists(self.fai_file):
            self.read()
        else:
            self.fai_file = f"{fasta_file_name}.gz.fai"
            if os.path.exists(self.fai_file):
                self.read()
            else:
                raise G2GFastaError("Unable to find fasta index file")

    def __getitem__(self, entry):
        try:
            if isinstance(entry, int):
                entry = tuple(self.records.keys())[entry]

            return self.records[entry]
        except KeyError:
            raise G2GFastaError(f"{entry} not in {self.fai_file}.")

    def __iter__(self):
        for keys in self.records:
            yield keys

    def __repr__(self):
        return f'{self.__class__}("{self.fai_file}")'

    def read(self):
        with open(self.fai_file) as index:
            for line in index:
                line = line.strip()
                idx, length, offset, line_len, line_len_bytes = line.split("\t")
                if idx in self.records:
                    raise ValueError(f"A duplicate entry was found: {idx}")
                else:
                    self.records[idx] = FAIEntry(
                        idx,
                        int(length),
                        int(offset),
                        int(line_len),
                        int(line_len_bytes),
                    )

    def get_pos(self, seq_id, start, end) -> tuple[int, int, int]:
        """
        Get the byte positions of the fasta file using the specified location

        Args:
            seq_id: The sequence identifier.
            start: The start position.
            end: The end position.

        Returns:
            A tuple of starting byte, ending byte, and byte length.
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

    def dump(self):
        print("FAI")
        print(f"Fasta File: {self.fasta_file}")
        print(f"FAI File: {self.fai_file}")
        for record in self.records:
            print(str(self.records[record]))


def extract(
    fasta_file: str | FastaFile,
    locations: g2g.Region | list[g2g.Region],
    output_file_name: str | None = None,
    reverse: bool | None = False,
    complement: bool | None = False,
    raw: bool | None = False,
    debug_level: int | None = 0,
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
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).

    Raises:
        G2GRegionError: When a Region is not correct.
    """
    start_time = time.time()
    logger = g2g.get_logger(debug_level)

    if isinstance(fasta_file, pysam.FastaFile):
        fasta = fasta_file
    else:
        fasta_input = g2g_utils.check_file(fasta_file)
        fasta = pysam.FastaFile(fasta_input)

    logger.warn(f"FASTA FILE: {fasta.filename}")

    if output_file_name:
        output_file_name = g2g_utils.check_file(output_file_name, "w")
        fasta_out = open(output_file_name, "w")
    else:
        fasta_out = sys.stdout

    try:
        if not isinstance(locations, list):
            locations = [locations]

        for location in locations:
            logger.debug(f"LOCATION: {location}")
            if not isinstance(location, g2g.Region):
                raise G2GRegionError(f"Found: {location}, which is not a Region")

            if location.seq_id not in fasta.references:
                continue

            sequence = fasta.fetch(location.seq_id, location.start, location.end)

            if location.strand == "+":
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
                fasta_out.write("\n")
            else:
                # internally storing 0 base, output 1 base
                start = location.start + 1
                end = location.end
                len_seq = len(sequence)

                # if someone requests something to long, just truncate
                # the length
                if len_seq < end:
                    end = len_seq + start - 1

                fasta_id = f">{location.seq_id}:{start}-{end}"
                logger.debug(f"Fasta ID: {fasta_id}")
                logger.debug(f"{sequence[:10]}...{sequence[-10:]}")

                if location.name:
                    fasta_id = f">{location.name} {location.seq_id}:{start}-{end}"
                    logger.debug(f"Location name used, Fasta ID: {fasta_id}")

                fasta_out.write(f"{fasta_id}\n")
                g2g_utils.write_sequence(sequence, fasta_out)

    except G2GRegionError as e:
        logger.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(e.msg.rstrip())
        raise e

    fasta_out.close()
    fmt_time = g2g_utils.format_time(start_time, time.time())
    logger.warn(f"Execution complete: {fmt_time}")


def extract_id(
    fasta_file: str | FastaFile,
    identifier: str,
    output_file_name: str | None = None,
    reverse: bool | None = False,
    complement: bool | None = False,
    raw: bool | None = False,
    debug_level: int | None = 0,
):
    """
    Extract the Fasta sequences.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        identifier: An identifier.
        output_file_name: Name of output file, None for sys.stdout.
        reverse: True to reverse the extracted sequence.
        complement: True to complement the extracted sequence.
        raw: True for just the sequence to be output, no fasta line.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    logger = g2g.get_logger(debug_level)

    start_time = time.time()

    if isinstance(fasta_file, pysam.FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = pysam.FastaFile(fasta_file)

    logger.warn(f"FASTA FILE: {fasta.filename}")

    if output_file_name:
        output_file_name = g2g_utils.check_file(output_file_name, "w")
        fasta_out = open(output_file_name, "w")
    else:
        fasta_out = sys.stdout

    try:
        sequence = fasta[identifier]

        if reverse and not complement:
            sequence = g2g_utils.reverse_sequence(sequence)
        elif not reverse and complement:
            sequence = g2g_utils.complement_sequence(sequence)
        elif reverse and complement:
            sequence = g2g_utils.reverse_complement_sequence(sequence)

        if raw:
            fasta_out.write(sequence)
        else:
            fasta_id = f">{identifier}\n"
            fasta_out.write(fasta_id)

            for line in g2g_utils.wrap_sequence(sequence):
                fasta_out.write(line.strip())
                fasta_out.write("\n")

    except G2GRegionError as e:
        logger.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(e.msg.rstrip())
        raise e

    fasta_out.close()
    fmt_time = g2g_utils.format_time(start_time, time.time())
    logger.warn(f"Execution complete: {fmt_time}")


def diff_files(file_name_1: str, file_name_2: str) -> None:
    """
    Print the number of differences in the 2 files.

    Args:
        file_name_1: The first Fasta file name.
        file_name_2: The second Fasta file name.
    """
    fasta_1 = pysam.FastaFile(file_name_1)
    fasta_2 = pysam.FastaFile(file_name_2)

    fa_equal = {}
    fa_diff = {}

    for seq_id1 in fasta_1.references:
        if seq_id1 in fasta_2.references:
            print(f"Comparing {seq_id1}")

            if fasta_1.fetch(seq_id1) == fasta_2.fetch(seq_id1):
                fa_equal[seq_id1] = seq_id1
            else:
                fa_diff[seq_id1] = seq_id1

    print(f"# EQUAL: {len(fa_equal)}")
    print(f"# DIFF: {len(fa_diff)}")

    for seq in fa_diff:
        print(seq)


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

        if s[0] == " ":
            continue
        elif s[0] == "-":
            print(f'Delete "{s[-1]}" from position {i}')
        elif s[0] == "+":
            print(f'Add "{s[-1]}" to position {i}')
    print()


def get_pos(fai: FAI, chrom: str, start: int, end: int) -> tuple[int, int, int]:
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
        fasta_file_name: str,
        output_file_name: str = None,
        length: int = 60,
        debug_level: int = 0,
):
    """
    Reformat a Fasta file to specified maximum line length.

    Args:
        fasta_file_name: The name of the Fasta file.
        output_file_name: Name of the output file or None for stdout.
        length: Maximum line length.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    fasta_file_name = g2g_utils.check_file(fasta_file_name)
    output = sys.stdout
    logger = g2g.get_logger(debug_level)

    if output_file_name:
        output_file = g2g_utils.check_file(output_file_name, "w")
        output = open(output_file, "w")

    with open(fasta_file_name, "r") as fasta_fd:
        new_sequence = StringIO()
        new_header = ">UNKNOWN"

        for line in fasta_fd:
            line = line.strip()

            if line[0] == ">":
                if len(new_sequence.getvalue()) > 0:
                    output.write(new_header)
                    output.write("\n")
                    g2g_utils.write_sequence(new_sequence.getvalue(), output, length)

                    logger.info(f"Reformatting {new_header}")

                new_sequence = StringIO()
                new_header = line
            else:
                new_sequence.write(line)

        logger.info(f"Reformatting {new_header}")
        output.write(new_header)
        output.write("\n")
        g2g_utils.write_sequence(new_sequence.getvalue(), output, length)

    output.close()


def fasta_extract_transcripts(
    fasta_file: str | FastaFile,
    database_file_name: str,
    output: str | None,
    raw: bool | None = False,
    debug_level: int | None = 0,
) -> None:
    """
    Extract the transcripts sequences from the fasta_file given the database.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        database_file_name: Name of the database file.
        output: Name of the output file or None for stdout.
        raw: True to omit Fasta sequence name and/or info.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start = time.time()
    logger = g2g.get_logger(debug_level)

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file_name = g2g_utils.check_file(database_file_name)
    fasta_out = sys.stdout

    if output:
        output = g2g_utils.check_file(output, "w")
        fasta_out = open(output, "w")

    logger.warn(f"FASTA FILE: {fasta.filename}")
    logger.warn(f"DATABASE FILE: {database_file_name}")
    logger.warn(f"OUTPUT FILE: {fasta_out.name}")

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file_name)

        for i, transcript in enumerate(transcripts):
            logger.debug(f"Transcript={transcript}")

            if transcript.seqid not in fasta.references:
                logger.debug("Skipping, transcript seqid not in fasta references")
                logger.debug(f"{transcript.seqid} not in {fasta.references}")
                continue

            new_sequence = StringIO()
            locations = []
            for ensembl_id, exon in transcript.exons.items():
                logger.debug(f"Exon ID={ensembl_id};{exon}")

                partial_seq = fasta.fetch(exon.seqid, exon.start - 1, exon.end)
                partial_seq_str = str(partial_seq)

                if transcript.strand == 1:
                    new_sequence.write(partial_seq_str)
                else:
                    partial_seq_str = str(
                        g2g_utils.reverse_complement_sequence(partial_seq)
                    )
                    new_sequence.write(partial_seq_str)

                logger.debug(
                    f"{exon.seqid}:{exon.start}-{exon.end} "
                    f"(Length: {len(partial_seq)})\n{partial_seq_str}"
                )
                locations.append(f"{exon.seqid}:{exon.start}-{exon.end}")

            if raw:
                fasta_out.write(new_sequence.getvalue())
            else:
                t_strand = "-" if transcript.strand == -1 else "+"
                loc = "|".join(locations)
                fasta_id = f">{transcript.ensembl_id} {t_strand}|{loc}\n"
                fasta_out.write(fasta_id)

                for line in g2g_utils.wrap_sequence(new_sequence.getvalue()):
                    fasta_out.write(line.strip())
                    fasta_out.write("\n")

    except G2GValueError as e:
        logger.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(e.msg.rstrip())
        raise e

    fmt_time = g2g_utils.format_time(start, time.time())
    logger.warn(f"Execution complete: {fmt_time}")


def fasta_extract_exons(
    fasta_file: str | FastaFile,
    database_file_name: str,
    output: str | None,
    raw: bool | None = False,
    debug_level: int | None = 0,
) -> None:
    """
    Extract the exons sequences from the fasta_file given the database.

    Args:
        fasta_file: Either the name of the Fasta file or a FastaFile object.
        database_file_name: Name of the database file.
        output: Name of the output file or None for stdout.
        raw: True to omit Fasta sequence name and/or info.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start = time.time()
    logger = g2g.get_logger(debug_level)

    if isinstance(fasta_file, FastaFile):
        fasta = fasta_file
    else:
        fasta_file = g2g_utils.check_file(fasta_file)
        fasta = FastaFile(fasta_file)

    database_file_name = g2g_utils.check_file(database_file_name)
    fasta_out = sys.stdout

    if output:
        output = g2g_utils.check_file(output, "w")
        fasta_out = open(output, "w")

    logger.warn(f"FASTA FILE: {fasta.filename}")
    logger.warn(f"DATABASE FILE: {database_file_name}")
    logger.warn(f"OUTPUT FILE: {fasta_out.name}")

    try:
        transcripts = gtf_db.get_transcripts_simple(database_file_name)
        for i, transcript in enumerate(transcripts):
            if transcript.seqid not in fasta.references:
                continue

            for ensembl_id, exon in transcript.exons.items():
                logger.debug(f"Exon={exon}")

                partial_seq = fasta.fetch(exon.seqid, exon.start - 1, exon.end)
                partial_seq_str = partial_seq

                if transcript.strand == -1:
                    rev_seq = g2g_utils.reverse_complement_sequence(partial_seq)
                    partial_seq_str = str(rev_seq)

                logger.debug(
                    f"{exon.seqid}:{exon.start}-{exon.end} "
                    f"(Length: {len(partial_seq)})\n{partial_seq_str}"
                )

                if raw:
                    fasta_out.write(partial_seq_str)
                else:
                    fasta_id = (
                        f">{exon.ensembl_id} " f"{exon.seqid}:{exon.start}-{exon.end}\n"
                    )
                    fasta_out.write(fasta_id)

                    for line in g2g_utils.wrap_sequence(partial_seq_str):
                        fasta_out.write(line.strip())
                        fasta_out.write("\n")

    except G2GValueError as e:
        logger.info(e.msg.rstrip())
        raise e
    except G2GFastaError as e:
        logger.info(e.msg.rstrip())
        raise e

    fmt_time = g2g_utils.format_time(start, time.time())
    logger.warn(f"Execution complete: {fmt_time}")

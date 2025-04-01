"""
Collection of functions related to BED files.

This module provides utilities for working with BED (Browser Extensible Data) files,
including parsing, iteration, and coordinate conversion. BED files use 0-based,
half-open coordinate system.
"""

# standard library imports
import collections
import sys
import time
from typing import Any, Iterator, TextIO

# 3rd party library imports
# none

# local library imports
from g2gtools.exceptions import G2GBedError, G2GError, KeyboardInterruptError
from g2gtools.vci import VCIFile
import g2gtools.g2g_utils as g2g_utils

# define BED record structure
bed_fields = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'extra']
BEDRecord = collections.namedtuple('BEDRecord', bed_fields)

# module logger
logger = g2g_utils.get_logger('g2gtools')


class BED:
    """
    Simple BED object for parsing and iterating BED files.

    This class provides functionality to read and iterate through BED files,
    supporting transparent gzip decompression. It parses each line of a BED file
    into a BEDRecord object with appropriate fields.

    Attributes:
        file_name (str): Path to the BED file.
        current_line (str | None): The current line being processed.
        current_line_is_bed (bool): Whether the current line is a BED record.
        current_record (BEDRecord | None): The current parsed BED record.
        reader: File reader object for the BED file.
        n_items (int | None): Number of fields in each BED record.
        current_line_no (int): Current line number being processed.
    """

    # Type annotations for instance variables
    file_name: str
    current_line: str | None
    current_line_is_bed: bool
    current_record: BEDRecord | None
    reader: TextIO  # File-like object
    n_items: int | None
    current_line_no: int

    def __init__(self, file_name: str) -> None:
        """
        Initialize a new BED file reader.

        Args:
            file_name: Path to the BED file to read.
        """
        self.file_name = file_name
        self.current_line = None
        self.current_line_is_bed = False
        self.current_record = None
        self.reader = g2g_utils.open_resource(file_name)
        self.n_items = None
        self.current_line_no = 0

    def __iter__(self) -> 'BED':
        """
        Make the BED object iterable.

        Returns:
            The BED object itself as an iterator.
        """
        return self

    def __next__(self) -> BEDRecord | None:
        """
        Get the next BED record from the file.

        Reads the next line from the BED file, parses it, and returns a BEDRecord
        object. Skips empty lines and handles track lines.

        Returns:
            A BEDRecord object representing the next record in the file,
            or None for track lines.

        Raises:
            G2GBedError: When there is an improperly formatted BED file.
            StopIteration: When the end of the file is reached.
        """
        try:
            self.current_line = g2g_utils.s(next(self.reader))
            self.current_line_no += 1

            # Skip empty lines
            while self.current_line and len(self.current_line.strip()) == 0:
                self.current_line = g2g_utils.s(next(self.reader))
                self.current_line_no += 1

            # Handle track lines
            if self.current_line.startswith('track'):
                self.current_line = self.current_line.strip()
                self.current_line_is_bed = False
                self.current_record = None
                return None

            self.current_line_is_bed = True
            elem = self.current_line.strip().split('\t')

            # Check consistent number of fields
            if not self.n_items:
                self.n_items = len(elem)
            elif self.n_items != len(elem):
                raise G2GBedError(
                    f'Inconsistent number of fields in BED file at line {self.current_line_no}'
                )

            try:
                bed_data = {
                    'chrom': elem[0],
                    'start': int(elem[1]),
                    'end': int(elem[2]),
                    'name': elem[3] if self.n_items > 3 else None,
                    'score': elem[4] if self.n_items > 4 else None,
                    'strand': elem[5] if self.n_items > 5 else None,
                    'extra': elem[6:] if self.n_items > 6 else None,
                }
                self.current_record = BEDRecord(**bed_data)
                return self.current_record
            except IndexError:
                raise G2GBedError(
                    f'Improperly formatted BED file, line number: {self.current_line_no}, '
                    f'line: {self.current_line}'
                )
            except ValueError:
                raise G2GBedError(
                    f'Improperly formatted BED file, line number: {self.current_line_no}, '
                    f'line: {self.current_line}'
                )
        except StopIteration:
            raise

    def close(self) -> None:
        """
        Close the BED file reader.

        This method should be called when done with the BED file to release resources.
        """
        if hasattr(self.reader, 'close'):
            self.reader.close()


def convert_bed_file(
        vci_file: str | VCIFile,
        bed_file_name_in: str,
        bed_file_name_out: str | None = None,
        reverse: bool = False,
) -> None:
    """
    Convert BED coordinates between genome assemblies using a VCI file.

    This function transforms coordinates in a BED file from one genome assembly
    to another using the mapping information in a VCI file. Unmapped regions
    are written to a separate file.

    Args:
        vci_file: Name of the VCI file or a VCIFile object.
        bed_file_name_in: Path to the input BED file to convert.
        bed_file_name_out: Path to the output BED file. If None, output is written to stdout.
        reverse: Whether to process VCI mappings in reverse direction.

    Raises:
        G2GError: If an error occurs during conversion.
        KeyboardInterruptError: If the user interrupts the process.
    """
    start_time = time.time()
    bed_out: TextIO | None = None
    bed_unmapped_file: TextIO | None = None

    try:
        # initialize VCI file if needed
        if not isinstance(vci_file, VCIFile):
            vci_path = g2g_utils.check_file(vci_file)
            vci_file = VCIFile(vci_path)
            logger.warning(f'Input VCI File: {vci_file.filename}')
            vci_file.parse(reverse)
        else:
            logger.warning(f'Input VCI File: {vci_file.filename}')

        # check input BED file
        bed_file_name_in = g2g_utils.check_file(bed_file_name_in)
        logger.warning(f'Input BED File: {bed_file_name_in}')

        # set up output files
        if bed_file_name_out:
            bed_file_name_out = g2g_utils.check_file(bed_file_name_out, 'w')
            unmapped_file = f'{bed_file_name_out}.unmapped'
            bed_out = open(bed_file_name_out, 'w')
            bed_unmapped_file = open(unmapped_file, 'w')
            logger.info(f'Output BED File: {bed_file_name_out}')
            logger.info(f'Output UNMAPPED File: {unmapped_file}')
        else:
            input_dir, input_name = g2g_utils.get_dir_and_file(bed_file_name_in)
            unmapped_file = f'{input_name}.unmapped'
            unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
            bed_out = sys.stdout
            bed_unmapped_file = open(unmapped_file, 'w')
            logger.info('Output BED File: stdout')
            logger.info(f'Output UNMAPPED File: {unmapped_file}')

        # log VCI file properties
        logger.info(f'VCI File is diploid: {vci_file.is_diploid()}')
        left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

        # process BED file
        logger.warning('Converting BED file...')
        bed_file = BED(bed_file_name_in)

        # counters for statistics
        total = 0
        success = 0
        fail = 0

        # process each BED record
        for record in bed_file:
            if record is None:  # Skip track lines
                continue

            logger.debug(f'ORIGINAL: {str(bed_file.current_line).strip()}')
            total += 1

            if total % 100000 == 0:
                logger.info(f'Processed {total:,} lines')

            for lr in left_right:
                seq_id = f'{record.chrom}{lr}'

                # find mappings for this region
                mappings = vci_file.find_mappings(
                    seq_id, record.start - 1, record.end
                )

                # handle unmapped regions
                if mappings is None:
                    logger.debug('\tFail due to no mappings')
                    bed_unmapped_file.write(bed_file.current_line)
                    fail += 1
                    continue

                logger.debug(f'{len(mappings)} mappings found')
                success += 1

                # calculate new coordinates
                start = mappings[0].to_start + 1
                end = mappings[-1].to_end

                # create new BED record
                elem = bed_file.current_line.rstrip().split('\t')
                logger.debug(f'({record.start - 1}, {record.end})=>({start}, {end})')
                logger.debug(elem)

                elem[0] = seq_id
                elem[1] = str(start)  # convert to string for joining
                elem[2] = str(end)  # convert to string for joining

                temp_elem = '\t'.join(map(str, elem))
                logger.debug(f'     NEW: {temp_elem}')
                bed_out.write(f'{temp_elem}\n')

        # close files
        if bed_file_name_out and bed_out:
            bed_out.close()

        if bed_unmapped_file:
            bed_unmapped_file.close()

        bed_file.close()

        # log statistics
        logger.warning(f'Converted {success:,} of {total:,} records')
        logger.warning('BED File converted')

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        # log execution time
        fmt_time = g2g_utils.format_time(start_time, time.time())
        logger.warning(f'Time: {fmt_time}')
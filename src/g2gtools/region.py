# standard library imports
import argparse
import re
import sys

# 3rd party library imports
# none

# local library imports
from g2gtools.exceptions import G2GRegionError

# from g2gtools.vci import VCIFile
# from g2gtools.g2g_utils import get_logger

# logger = get_logger('g2gtools')

REGEX_REGION = re.compile(
    r'([\w.-]*)\s*(?:-|:)?\s*(\d+)\s*(MB|M|K)?\s*(?:-|:)?\s*(\d+)?\s*(MB|M|K)?',
    re.IGNORECASE,
)

#REGEX_REGION_CHR = re.compile(
#    r'(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?',
#    re.IGNORECASE,
#)


def exit(message: str = '', parser: argparse.ArgumentParser = None) -> None:
    """
    Print message, help, and exit.

    Args:
        message: The message to print.
        parser: The argument parser.
    """
    if parser:
        if message:
            parser.error(message)
        else:
            parser.error('Unknown error specified')
    else:
        if message:
            sys.stderr.write(f'[g2gtools] Error: {message}\n')

    sys.exit(1)


#
#
# Region
#
#

class Region:
    """
    Represents a genomic region with sequence ID, coordinates, and strand information.

    This class encapsulates a genomic region defined by a sequence identifier,
    start and end positions, and strand orientation. It provides validation for
    coordinates and strand values, and supports both 0-based and 1-based coordinate
    systems.

    Attributes:
        seq_id (str): Sequence identifier (e.g., chromosome name).
        _strand (str | None): Strand orientation ('+' or '-').
        _start (int | None): Start position of the region.
        _end (int | None): End position of the region.
        original_base (int): Coordinate system base (0 or 1).
        name (str | None): Optional name for the region.
    """
    # Type annotations for instance variables
    seq_id: str
    _strand: str | None
    _start: int | None
    _end: int | None
    original_base: int
    name: str | None

    def __init__(
            self,
            seq_id: str,
            start: int | None = None,
            end: int | None = None,
            strand: str | int = '+',
            name: str | None = None,
            original_base: int = 0,
    ) -> None:
        """
        Initialize a new Region instance.

        Args:
            seq_id: Sequence identifier (e.g., chromosome name).
            start: Start position of the region. Defaults to None.
            end: End position of the region. If None and start is provided,
                 defaults to start + 1. Defaults to None.
            strand: Strand orientation ('+', '-', 1, or -1). Defaults to '+'.
            name: Optional name for the region. Defaults to None.
            original_base: Coordinate system base (0 or 1). Defaults to 0.
        """
        if start:
            end = end if end else start + 1
        self.seq_id = seq_id
        self._strand = None
        self._start = None
        self._end = None
        self.original_base = original_base
        self.strand = strand
        self.start = start
        self.end = end
        self.name = name

    @property
    def strand(self) -> str | None:
        """
        Get the strand orientation.

        Returns:
            The strand orientation ('+' or '-') or None if not set.
        """
        return self._strand

    @strand.setter
    def strand(self, value: str | int | None) -> None:
        """
        Set the strand orientation.

        Args:
            value: Strand value ('+', '-', 1, -1, or None).

        Raises:
            G2GRegionError: If the strand value is invalid.
        """
        if value is not None:
            if value in ('+', '1', 1):
                self._strand = '+'
            elif value in ('-', '-1', -1):
                self._strand = '-'
            else:
                raise G2GRegionError(
                    f'Illegal value for strand {value}, must be +, -, 1, -1'
                )

    @property
    def start(self) -> int | None:
        """
        Get the start position.

        Returns:
            The start position or None if not set.
        """
        return self._start

    @start.setter
    def start(self, value: int | None) -> None:
        """
        Set the start position.

        Args:
            value: Start position value or None.

        Raises:
            G2GRegionError: If the start position is negative.
        """
        if value is not None:
            if value < 0:
                raise G2GRegionError(
                    f'Illegal value for start {value}, start must be >= 0'
                )
            self._start = value

    def get_start(self) -> int | None:
        """
        Get the start position adjusted for the coordinate system.

        Returns:
            The start position adjusted for 0-based or 1-based coordinates,
            or None if start is not set.
        """
        if self.start is None:
            return None
        if self.original_base == 1:
            return self.start
        return self.start - 1

    @property
    def end(self) -> int | None:
        """
        Get the end position.

        Returns:
            The end position or None if not set.
        """
        return self._end

    @end.setter
    def end(self, value: int | None) -> None:
        """
        Set the end position.

        If value is None and start is set, end will be set to start + 1.

        Args:
            value: End position value or None.

        Raises:
            G2GRegionError: If the end position is less than or equal to the start position.
        """
        if value is not None:
            if self.start is not None and value <= self.start:
                raise G2GRegionError(
                    f'Illegal value for end {value}, end must be >= start {self.start}'
                )
            self._end = value
        else:
            if self.start is not None:
                self._end = self.start + 1

    def __str__(self) -> str:
        """
        Get a string representation of the region.

        Returns:
            A string representation including name (if available) and coordinates.
        """
        # Adjust start position for display based on coordinate system
        display_start = self.start + 1 if self.original_base == 0 else self.start

        if self.name:
            return f'{self.name} {self.seq_id}:{display_start}-{self.end}'
        else:
            return f'{self.seq_id}:{display_start}-{self.end}'

    def __repr__(self) -> str:
        """
        Get a representation of the region for debugging.

        Returns:
            A string representation including coordinates and strand.
        """
        display_start = self.start + 1 if self.original_base == 0 else self.start
        return f'{self.seq_id}:{display_start}-{self.end} ({self.strand})'

def get_multiplier(factor: str) -> int:
    """
    Convert a size factor string into its numeric multiplier.

    Args:
        factor: The string 'mb', 'm', or 'k' (case-insensitive)

    Returns:
        The corresponding multiplier: 1000000 for 'mb'/'m', 1000 for 'k',
        or 1 otherwise
    """
    factor = factor.lower() if factor else ''
    if factor == 'mb':
        return 1000000  # Fixed: MB is 1,000,000 not 10,000,000
    elif factor == 'm':
        return 1000000
    elif factor == 'k':
        return 1000
    return 1


def parse_region(
        location_str: str,
        base: int = 0,
        name: str | None = None
) -> Region:
    """
    Parse a string and return a Region object.

    Format expected is seqid:start-end or seqid:start or similar variations,
    where ':' and '-' can be other delimiters. Supports 'K', 'M', and 'MB' suffixes.

    Args:
        location_str: A string representing a genomic location.
        base: Coordinate system base (0 or 1). Defaults to 0.
        name: Optional name for the region. Defaults to None.

    Returns:
        A Region object representing the parsed location.

    Raises:
        G2GRegionError: If the location string cannot be parsed.
    """
    # clean up input
    clean_location = str(location_str).strip().replace(',', '')

    # match the location pattern
    loc_match = REGEX_REGION.match(clean_location)
    if not loc_match:
        raise G2GRegionError(f"Cannot parse location '{location_str}'")

    # extract components
    loc_groups = loc_match.groups()
    identifier = loc_groups[0]
    start_base_str = loc_groups[1]
    start_mult = loc_groups[2]
    end_base_str = loc_groups[3]
    end_mult = loc_groups[4]

    # validate identifier
    if not identifier:
        raise G2GRegionError('Missing sequence identifier')

    # parse start position
    if not start_base_str:
        raise G2GRegionError('Missing start position')

    try:
        start_base = int(start_base_str)
    except ValueError:
        raise G2GRegionError('Start position is not numeric')

    # apply start multiplier
    start = start_base
    if start_mult:
        start_mult = start_mult.lower()
        if start_mult not in ['mb', 'm', 'k']:
            raise G2GRegionError(f'Unknown quantifier: {start_mult}')
        start = start_base * get_multiplier(start_mult)
    elif base == 1:
        # Adjust for 1-based coordinates
        start -= 1

    # parse end position
    end = None
    if end_base_str:
        try:
            end_base = int(end_base_str)
        except ValueError:
            raise G2GRegionError('End position is not numeric')

        # apply end multiplier
        if end_mult:
            end_mult = end_mult.lower()
            if end_mult not in ['mb', 'm', 'k']:
                raise G2GRegionError(f'Unknown quantifier: {end_mult}')
            end = end_base * get_multiplier(end_mult)
        else:
            end = end_base

    # create and return the Region
    return Region(identifier, start, end, '+', original_base=base, name=name)

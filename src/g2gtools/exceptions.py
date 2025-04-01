"""
Collection of module errors.

This module defines the exception hierarchy used throughout the g2gtools package.
All custom exceptions inherit from the base G2GError class, which itself inherits
from the standard Exception class.
"""


class G2GError(Exception):
    """
    Base exception class for all g2gtools errors.

    This is the root exception from which all other g2gtools-specific
    exceptions inherit.
    """
    def __init__(self, msg: str | None = None) -> None:
        """
        Initialize a new G2GError.

        Args:
            msg: The error message. Defaults to None.
        """
        # no need to store msg as an attribute since Exception already does this
        super().__init__(msg)


# File-related exceptions
class G2GBAMError(G2GError):
    """Exception raised for errors related to BAM/SAM files."""
    pass


class G2GBedError(G2GError):
    """Exception raised for errors related to BED files."""
    pass


class G2GFastaError(G2GError):
    """Exception raised for errors related to FASTA sequence files."""
    pass


class G2GVCFError(G2GError):
    """Exception raised for errors related to Variant Call Format files."""
    pass


class G2GGTFError(G2GError):
    """Exception raised for errors related to Gene Transfer Format files."""
    pass


# Format and parsing exceptions
class G2GCigarFormatError(G2GError):
    """Exception raised for errors in CIGAR string format in SAM/BAM files."""
    pass


class G2GValueError(G2GError):
    """
    Exception raised for errors in the value of parameters.

    This exception is used when a function receives a parameter with an
    inappropriate value (e.g., incorrect data type, out of range).
    """
    pass


# Operation-specific exceptions
class G2GFetchError(G2GError):
    """Exception raised for errors during sequence or data fetching operations."""
    pass


class G2GRegionError(G2GError):
    """Exception raised for errors related to genomic region specifications or operations."""
    pass


class KeyboardInterruptError(Exception):
    """
    Exception raised to handle keyboard interrupts.

    This exception is used to properly catch and handle keyboard
    interrupts (Ctrl+C) during long-running operations.
    """
    pass

# -*- coding: utf-8 -*-

#
# Collection of module errors
#


class G2GError(Exception):
    """
    Simple base exception, root of all G2G exceptions
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GError, self).__init__(self.msg)


class G2GValueError(G2GError):
    """
    Generic value error for G2G
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GValueError, self).__init__(self.msg)


class G2GBAMError(G2GError):
    """
    BAM/SAM file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GBAMError, self).__init__(self.msg)


class G2GBedError(G2GError):
    """
    BED file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GBedError, self).__init__(self.msg)


class G2GChainFileError(G2GError):
    """
    Chain file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GChainFileError, self).__init__(self.msg)


class G2GCigarFormatError(G2GError):
    """
    Cigar exception/error
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GCigarFormatError, self).__init__(self.msg)


class G2GFastaError(G2GError):
    """
    Fasta file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GFastaError, self).__init__(self.msg)


class G2GFetchError(G2GError):
    """
    G2G Fetch errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GFetchError, self).__init__(self.msg)


class G2GLocationError(G2GError):
    """
    Location errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GLocationError, self).__init__(self.msg)


class G2GVCFError(G2GError):
    """
    VCF file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GVCFError, self).__init__(self.msg)


class G2GGTFError(G2GError):
    """
    GTF file errors
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(G2GGTFError, self).__init__(self.msg)



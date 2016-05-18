# -*- coding: utf-8 -*-

from .exceptions import G2GError, G2GBAMError, G2GFastaError, G2GFetchError, G2GCigarFormatError, G2GChainFileError,\
    G2GBedError, G2GVCFError, G2GLocationError, G2GValueError, G2GGTFError
from .g2g import file_convert, offset2chain, gtf2chain, db2chain, fasta_transform
from .g2g_commands import command_convert, command_offset2chain, command_gtf2chain, command_gtf2db, command_db2chain, \
    command_vcf2chain, command_fasta_extract, command_fasta_patch, command_fasta_transform, command_parse_location
from .g2g_utils import parse_location, complement_sequence, reverse_sequence, reverse_complement_sequence
from .gtf_db import gtf2db
from .chain import ChainFile
from .fasta import extract as fasta_extract
from .fasta_patch import fasta_patch
from .vcf2chain import vcf2chain
from .bins import bins


__version__ = '0.1.29'
__author__ = 'Matthew Vincent, The Jackson Laboratory'
__email__ = 'matt.vincent@jax.org'


__all__ = ['G2GError', 'G2GBAMError', 'G2GFastaError', 'G2GFetchError', 'G2GCigarFormatError', 'G2GChainFileError',
           'G2GBedError', 'G2GVCFError', 'G2GLocationError', 'G2GValueError', 'G2GGTFError',
           'file_convert', 'offset2chain', 'gtf2chain', 'gtf2db', 'db2chain', 'fasta_transform',
           'command_convert', 'command_offset2chain', 'command_gtf2chain', 'command_gtf2db',
           'command_db2chain', 'command_vcf2chain', 'command_fasta_extract', 'command_fasta_transform',
           'command_fasta_patch', 'command_parse_location',
           'parse_location', 'complement_sequence', 'reverse_sequence', 'reverse_complement_sequence',
           'gtf2db',
           'ChainFile',
           'fasta_extract',
           'fasta_patch',
           'vcf2chain'
           'bins']


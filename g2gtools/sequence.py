# -*- coding: utf-8 -*-

#
# Collection of functions related to Sequences
#

import re
from string import maketrans

BASES = re.compile(r"([ACTGNactgn]+)")
TRANS = maketrans('ACTGNactgn', 'TGACNtgacn')


def reverse_sequence(sequence):
    """
    Return the reverse of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse string or Sequence object
    """
    if isinstance(sequence, Sequence):
        return Sequence(seq=sequence.seq[::-1], id=sequence.id, start=sequence.start, end=sequence.end, strand=sequence.strand)
    elif isinstance(sequence, str):
        return sequence[::-1]
    else:
        raise ValueError("Cannot complement object of {0}, expecting string or Sequence object".format(type(sequence)))


def complement_sequence(sequence):
    """
    Return the complement of sequence

    :param sequence: either a string or Sequence object
    :return: the complement string or Sequence object
    """
    val = str(sequence)

    if len(re.findall(BASES, val)) == 1:  # finds invalid characters if > 1
        val = str(sequence).translate(TRANS)
    else:
        matches = re.findall(BASES, val)
        position = len(matches[0])
        raise ValueError("Sequence contains non-DNA character '{0}' at position {1:n}\n".format(val[position], position + 1))

    if isinstance(sequence, Sequence):
        return Sequence(seq=val, id=sequence.id, start=sequence.start, end=sequence.end, strand=sequence.strand)
    elif isinstance(sequence, str):
        return val
    else:
        raise ValueError("Cannot complement object of {0}, expecting string or Sequence object".format(type(sequence)))


def reverse_complement_sequence(sequence):
    """
    Return the reverse-complement of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse-complement string or Sequence object
    """
    return reverse_sequence(complement_sequence(sequence))


class Sequence(object):
    """
    Simple class to abstract sequence information.

    Simplest form is:
        seq = Sequence("ACGTACGT")

    Can also take the following keyword arguments:

    id = identifier
    seq = string of ACGT
    start, end = positions
    strand
    """
    def __init__(self, *args, **kwargs):
        self.seq = None
        self.id = None
        self.start = None
        self.end = None
        self.strand = None

        #args -- tuple of anonymous arguments
        #kwargs -- dictionary of named arguments

        if 'seq' in kwargs and len(args) == 1:
            raise ValueError("Incorrect values passed to {0}".format(self.__class__))

        if len(args) == 1:
            self.seq = args[0]

        if 'seq' in kwargs:
            self.seq = kwargs['seq']
        if 'id' in kwargs:
            self.id = kwargs['id']
        if 'start' in kwargs:
            self.start = kwargs['start']
        if 'end' in kwargs:
            self.end = kwargs['end']
        if 'strand' in kwargs:
            self.strand = kwargs['strand']

    def __getitem__(self, n):
        if isinstance(n, slice):
            start, stop, step = n.indices(len(self))
            if stop == -1:
                stop = start
            else:
                stop = len(self) - stop
            if self.start and self.end:
                return self.__class__(seq=self.seq[n.start:n.stop:n.step], id=self.id,
                                  start=self.start + start, end=self.end - stop, strand=self.strand)
            else:
                return self.__class__(seq=self.seq[n.start:n.stop:n.step], id=self.id, strand=self.strand)
        elif isinstance(n, int):
            if n < 0:
                n += len(self)
            if self.start:
                return self.__class__(seq=self.seq[n], id=self.id, start=self.start + n, end=self.start + n, strand=self.strand)
            else:
                return self.__class__(seq=self.seq[n], id=self.id, strand=self.strand)

    def __str__(self):
        return self.seq

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __len__(self):
        return len(self.seq)

    @property
    def info(self):
        """ Return the fancy name for the sequence, including start, end
        """
        if len(self) > 50:
            # Shows the last three letters as it is often useful to see if there
            # is a stop codon at the end of a sequence.
            # Note total length is 44+3+3=50
            return "%s(id='%s', '%s...%s')" % (self.__class__.__name__, self.id, str(self)[:44], str(self)[-3:])

        return "%s(id='%s', '%s')" % (self.__class__.__name__, self.id, repr(self.seq))



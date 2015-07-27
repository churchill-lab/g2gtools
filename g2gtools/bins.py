# -*- coding: utf-8 -*-

# http://genomewiki.ucsc.edu/index.php/Bin_indexing_system

BIN_NEXT_SHIFT = 3    # How much to shift to get to finest bin
BIN_FIRST_SHIFT = 17  # How much to shift to get to next larger bin
#BIN_OFFSETS = [32768+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0]
BIN_OFFSETS = [4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0]
#BIN_OFFSETS = [512+64+8+1, 64+8+1, 8+1, 1, 0]
# for BED (0-based, half-open) or GFF (1-based, closed intervals)
COORD_OFFSETS = {'bed': 0, 'gff': 1}

def bins(start, stop, fmt='gff', one=True):
    """
    Uses the definition of a "genomic bin" described in Fig 7 of
    http://genome.cshlp.org/content/12/6/996.abstract.

    Parameters
    ----------
    one : boolean
        If `one=True` (default), then only return the smallest bin that
        completely contains these coordinates (useful for assigning a single
        bin).

        If `one=False`, then return the set of *all* bins that overlap these
        coordinates (useful for looking for features that could intersect)

    fmt : 'gff' | 'bed'
        This specifies 1-based start coords (gff) or 0-based start coords (bed)
    """
    # Jump to highest resolution bin that will fit these coords (depending on
    # whether we have a BED or GFF-style coordinate).
    start = (start - COORD_OFFSETS[fmt]) >> BIN_FIRST_SHIFT
    stop = (stop) >> BIN_FIRST_SHIFT

    # We always at least fit within the chrom, which is bin 1.
    bins = set([1])

    for offset in BIN_OFFSETS:
        # Since we're going from smallest to largest bins, the first one where
        # the feature's start and stop positions are both within the same bin
        # is the smallest one these coords fit within.
        if one:
            if start == stop:
                # Note that at this point, because of the bit-shifting, `start`
                # is the number of bins (at this current level).  So we need to
                # add it to `offset` to get the actual bin ID.
                return offset + start

        # See the Fig 7 reproduction above to see why range().
        bins.update(list(range(offset + start, offset + stop + 1)))

        # Move to the next level (8x larger bin size; i.e., 2**NEXT_SHIFT
        # larger bin size)
        start >>= BIN_NEXT_SHIFT
        stop >>= BIN_NEXT_SHIFT
    return bins


def bin_info():
    """
    Useful for debugging: how large is each bin, and what are the bin IDs?
    """
    print "level\t#bins\tstart\tend\tsize"
    level = -1
    for i, offset in reversed(list(enumerate(BIN_OFFSETS))):
        level += 1
        number_of_bins = 8 ** level

        binstart = offset
        try:
            binstop = BIN_OFFSETS[i+1]
            binstop = number_of_bins + binstart - 1
        except IndexError:
            binstop = binstart

        bin_size = 2 ** (BIN_FIRST_SHIFT + (i * BIN_NEXT_SHIFT))
        actual_size = bin_size

        # nice formatting
        bin_size, suffix = bin_size / 1024, 'Kb'
        if bin_size >= 1024:
            bin_size, suffix = bin_size / 1024, 'Mb'
        if bin_size >= 1024:
            bin_size, suffix = bin_size / 1024, 'Gb'
        size = '(%s %s)' % (bin_size, suffix)
        actual_size = '%s bp' % (actual_size)

        print('{level}\t{number_of_bins}\t{binstart}\t{binstop}\t{actual_size}\t{size}'.format(**locals()))


if __name__ == "__main__":
    bin_info()


#
# Collection of functions related to BAM and SAM files
#
# pysam uses 0-based coordinates

# standard library imports
from collections import namedtuple
import os
import pathlib
import re

# 3rd party library imports
import pysam

# local library imports
from . import __version__ as g2g_version
from . import g2g
from . import g2g_utils
from . import vci
from .exceptions import G2GBAMError
from .exceptions import G2GCigarFormatError


FLAG_NONE = 0x0             # base value
FLAG_PAIRED = 0x1           # template having multiple segments in sequencing
FLAG_PROPER_PAIR = 0x2      # each segment properly aligned according to the aligner
FLAG_UNMAP = 0x4            # segment unmapped
FLAG_MUNMAP = 0x8           # next segment in the template unmapped (mate unmapped)
FLAG_REVERSE = 0x10         # SEQ being reverse complemented
FLAG_MREVERSE = 0x20        # SEQ of the next segment in the template being reversed
FLAG_READ1 = 0x40           # the first segment in the template
FLAG_READ2 = 0x80           # the last segment in the template
FLAG_SECONDARY = 0x100      # secondary alignment
FLAG_QCFAIL = 0x200         # not passing quality controls
FLAG_DUP = 0x400            # PCR or optical duplicate
FLAG_SUPPLEMENTARY = 0x800  # supplementary alignment

REGEX_CIGAR = re.compile("(\d+)([\w=])")
REGEX_CIGAR_LENGTH = re.compile("\D")

CIGAR_M = 'M'
CIGAR_I = 'I'
CIGAR_D = 'D'
CIGAR_N = 'N'
CIGAR_S = 'S'
CIGAR_H = 'H'
CIGAR_P = 'P'
CIGAR_E = '='
CIGAR_X = 'X'

CIGAR_m = 0
CIGAR_i = 1
CIGAR_d = 2
CIGAR_n = 3
CIGAR_s = 4
CIGAR_h = 5
CIGAR_p = 6
CIGAR_e = 7
CIGAR_x = 8

CIGAR_N2C = {
    0: 'M',  # alignment match (can be a sequence match or mismatch)
    1: 'I',  # insertion to the reference
    2: 'D',  # deletion from the reference
    3: 'N',  # skipped region from the reference
    4: 'S',  # soft clipping (clipped sequences present in SEQ)
    5: 'H',  # hard clipping (clipped sequences NOT present in SEQ)
    6: 'P',  # padding (silent deletion from padded reference)
    7: '=',  # sequence match
    8: 'X',  # sequence mismatch
    '0': 'M',
    '1': 'I',
    '2': 'D',
    '3': 'N',
    '4': 'S',
    '5': 'H',
    '6': 'P',
    '7': '=',
    '8': 'X'
}

CIGAR_C2N = {
    'M': 0,
    'I': 1,
    'D': 2,
    'N': 3,
    'S': 4,
    'H': 5,
    'P': 6,
    '=': 7,
    'X': 8
}

Cigar = namedtuple("Cigar", ["code", "length", "start", "end"])


def open_alignment_file(file_name: str) -> pysam.AlignmentFile:
    """
    Open a pysam.AlignmentFile based upon the type of file it is.

    Args:
        file_name: Name of the file to open.

    Returns:
        An AlignmentFile.
    """
    # function to return the file extension
    file_extension = pathlib.Path(file_name).suffix
    print("File Extension: ", file_extension)

    alignment_file = None

    if file_extension.lower() == ".bam":
        alignment_file = pysam.AlignmentFile(file_name, "rb")
    elif file_extension.lower() == ".sam":
        alignment_file = pysam.AlignmentFile(file_name, "r")
    elif file_extension.lower() == ".cram":
        alignment_file = pysam.AlignmentFile(file_name, "rc")

    return alignment_file


def convert_bcsam_file(
        vci_file: str | vci.VCIFile,
        bcsam_file_name_in: str,
        bcsam_file_name_out: str,
        reverse: bool | None = False,
        debug_level: int | None = 0
) -> None:
    """
    Convert genome coordinates (in BAM/CRAM/SAM format) between assemblies.

    Args
        vci_file: Name of the VCI file or a VCIFile object.
        bcsam_file_name_in: Input BAM/CRAM/SAM file to convert.
        bcsam_file_name_out: Name of output BAM/CRAM/SAM file.
        reverse: True to process VCI in reverse.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    logger = g2g.get_logger(debug_level)

    if bcsam_file_name_out is None:
        raise G2GBAMError("Conversion of BAM/CRAM/SAM file needs output file")

    if not isinstance(vci_file, vci.VCIFile):
        vci_file = g2g_utils.check_file(vci_file)

    bcsam_file_name_in = g2g_utils.check_file(bcsam_file_name_in)
    bcsam_file_name_out = g2g_utils.check_file(bcsam_file_name_out, "w")
    unmapped_file_name = f"{bcsam_file_name_out}.unmapped"

    logger.warn(f"VCI FILE: {vci_file}")
    logger.warn(f"INPUT BAM/CRAM/SAM FILE: {bcsam_file_name_in}")
    logger.warn(f"OUTPUT BAM/CRAM/SAM FILE: {bcsam_file_name_out}")
    logger.warn(f"UNMAPPED FILE: {unmapped_file_name}")

    if not isinstance(vci_file, vci.VCIFile):
        logger.info("Parsing vci file...")
        vci_file = vci.VCIFile(vci_file, reverse=reverse)
        vci_file.parse(reverse=reverse)
        logger.info("VCI file parsed")

    alignment_file_in = open_alignment_file(bcsam_file_name_in)

    logger.info("Converting BAM file")

    new_header = alignment_file_in.header.to_dict()

    # replace 'HD'
    new_header["HD"] = {"VN": 1.0, "SO": "coordinate"}

    # replace SQ
    tmp = []
    name_to_id = {}

    for ref_name in vci_file.contigs:
        tmp.append({"LN": vci_file.contigs[ref_name], "SN": ref_name})
        name_to_id[ref_name] = alignment_file_in.get_tid(ref_name)

    new_header["SQ"] = tmp

    if "PG" not in new_header:
        new_header["PG"] = []

    new_header["PG"].append({"ID": "g2gtools", "VN": g2g_version})

    if "CO" not in new_header:
        new_header["CO"] = []

    new_header["CO"].append(f"Original file: {bcsam_file_name_in}")
    new_header["CO"].append(f"VCI File: {vci_file.filename}")

    dir, temp_file_name = os.path.split(bcsam_file_name_out)
    parts = temp_file_name.split(".")
    ext = parts[-1]

    if ext.lower() == "bam":
        new_file = pysam.AlignmentFile(bcsam_file_name_out, "wb", header=new_header)
        new_file_unmapped = pysam.AlignmentFile(unmapped_file_name, "wb", template=alignment_file_in)
    elif ext.lower() == "sam":
        new_file = pysam.Samfile(bcsam_file_name_out, "wh", header=new_header)
        new_file_unmapped = pysam.AlignmentFile(unmapped_file_name, "wh", template=alignment_file_in)
    else:
        raise G2GBAMError("Unable to create new file based upon file extension")

    total = 0
    total_unmapped = 0
    total_fail_qc = 0

    map_statistics = {
        "total": 0,
        "fail_cannot_map": 0,
        "success_simple": 0,
        "success_complex": 0
    }

    map_statistics_pair = {
        "total": 0,
        "fail_cannot_map": 0,
        "success_1_fail_2_simple": 0,
        "success_1_fail_2_complex": 0,
        "success_1_simple_2_fail": 0,
        "success_1_simple_2_simple": 0,
        "success_1_simple_2_complex": 0,
        "success_1_complex_2_fail": 0,
        "success_1_complex_2_simple": 0,
        "success_1_complex_2_complex": 0
    }

    try:
        while True:
            if total and total % 10000 == 0:
                status_success = 0
                status_failed = 0

                for k, v in map_statistics_pair.items():
                    if k.startswith("success"):
                        status_success += v
                    elif k.startswith("fail"):
                        status_failed += v

                logger.info(f"Processed {total:,} reads, {status_success:,} successful, {status_failed:,} failed")

            alignment = next(alignment_file_in)

            alignment_new = pysam.AlignedRead()
            read_chr = alignment_file_in.getrname(alignment.tid)

            # READ ONLY

            # aend                  aligned reference position of the read on the reference genome
            # alen                  aligned length of the read on the reference genome.
            # positions             a list of reference positions that this read aligns to
            # qend                  end index of the aligned query portion of the sequence (0-based, exclusive)
            # qlen                  Length of the aligned query sequence
            # qqual                 aligned query sequence quality values
            # qstart                start index of the aligned query portion of the sequence (0-based, inclusive)
            # query                 aligned portion of the read and excludes any flanking bases that were soft clipped
            # rlen                  length of the read

            # TRUE / FALSE (setting effects flag)

            # is_paired             true if read is paired in sequencing
            # is_proper_pair        true if read is mapped in a proper pair
            # is_qcfail             true if QC failure
            # is_read1              true if this is read1
            # is_read2              true if this is read2
            # is_reverse            true if read is mapped to reverse strand
            # is_secondary          true if not primary alignment
            # is_unmapped           true if read itself is unmapped
            # mate_is_reverse       true is read is mapped to reverse strand
            # mate_is_unmapped      true if the mate is unmapped

            # SET

            # cigar                 cigar as list of tuples
            # cigarstring           alignment as a string
            # flag                  properties flag
            # mapq                  mapping quality
            # pnext                 the position of the mate
            # pos                   0-based leftmost coordinate
            # pnext                 the position of the mate
            # qname                 the query name
            # rnext                 the reference id of the mate
            # seq                   read sequence bases, including soft clipped bases
            # tid                   target id, contains the index of the reference sequence in the sequence dictionary

            # DON"T NEED TO SET or SHOULD WE SET?

            # qual                  read sequence base qualities, including soft clipped bases
            # tags                  the tags in the AUX field

            # NOTE: Some tool require TLEN to be set.
            # `samtools fixmate` will correct TLEN if needed.
            # tlen = (read2_start + cigar bases) - read1_start (positive for R1 and negative for R2)
            # tlen                  insert size

            total += 1

            logger.debug("~"*80)
            logger.debug("Converting {0} {1} {2} {3}".format(alignment.qname, read_chr, alignment.pos, alignment.cigarstring))

            if alignment.is_qcfail:
                logger.debug("\tFail due to qc of old alignment")
                new_file_unmapped.write(alignment)
                total_fail_qc += 1
                continue

            if alignment.is_unmapped:
                logger.debug("\tFail due to unmapped old alignment")
                new_file_unmapped.write(alignment)
                total_unmapped += 1
                continue

            if not alignment.is_paired:
                logger.debug("SINGLE END ALIGNMENT")
                map_statistics["total"] += 1

                alignment_new.seq = alignment.seq
                alignment_new.flag = FLAG_NONE
                alignment_new.mapq = alignment.mapq
                alignment_new.qname = alignment.qname
                alignment_new.qual = alignment.qual
                alignment_new.tags = alignment.tags

                read_start = alignment.pos
                read_end = alignment.aend
                read_strand = "-" if alignment.is_reverse else "+"

                mappings = vci_file.find_mappings(read_chr, read_start, read_end)

                # unmapped
                if mappings is None:
                    logger.debug("\tFail due to no mappings")
                    new_file_unmapped.write(alignment)
                    map_statistics["fail_cannot_map"] += 1

                elif len(mappings) == 1:
                    if alignment.is_reverse:
                        alignment_new.flag |= FLAG_REVERSE

                    alignment_new.tid = name_to_id[mappings[0].to_chr]
                    alignment_new.pos = mappings[0].to_start
                    alignment_new.cigar = alignment.cigar
                    new_file.write(alignment_new)

                    logger.debug(f"\tSuccess (simple): {alignment_new.pos} {alignment_new.cigarstring}")
                    map_statistics["success_simple"] += 1

                else:
                    logger.debug(f"MAPPINGS: {len(mappings)}")
                    for m in mappings:
                        logger.debug(f"> {m}")

                    if alignment.is_reverse:
                        alignment_new.flag |= FLAG_REVERSE

                    alignment_new.tid = name_to_id[mappings[0].to_chr]
                    alignment_new.pos = mappings[0].to_start
                    alignment_new.cigar = convert_cigar(alignment.cigar, read_chr, vci_file, alignment.seq, read_strand, alignment.pos)
                    new_file.write(alignment_new)

                    logger.debug(f"\tSuccess (complex): {alignment_new.pos} {alignment_new.cigarstring}")
                    map_statistics["success_complex"] += 1

            else:
                logger.debug("PAIRED END ALIGNMENT")
                map_statistics_pair["total"] += 1

                alignment_new.seq = alignment.seq
                alignment_new.flag = FLAG_PAIRED
                alignment_new.mapq = alignment.mapq
                alignment_new.qname = alignment.qname
                alignment_new.qual = alignment.qual
                alignment_new.tags = alignment.tags

                if alignment.is_read1:
                    alignment_new.flag |= FLAG_READ1
                if alignment.is_read2:
                    alignment_new.flag |= FLAG_READ2

                if alignment.is_reverse:
                    alignment_new.flag |= FLAG_REVERSE
                if alignment.mate_is_reverse:
                    alignment_new.flag |= FLAG_MREVERSE
                if alignment.is_proper_pair:
                    alignment_new.flag |= FLAG_PROPER_PAIR

                read1_chr = alignment_file_in.getrname(alignment.tid)
                read1_start = alignment.pos
                read1_end = alignment.aend
                read1_strand = "-" if alignment.is_reverse else "+"
                read1_mappings = vci_file.find_mappings(
                    read1_chr, read1_start, read1_end
                )

                if alignment.mate_is_unmapped:
                    alignment_new.flag |= FLAG_MUNMAP
                else:
                    read2_chr = alignment_file_in.getrname(alignment.rnext)
                    read2_start = alignment.pnext
                    read2_end = read2_start + 1
                    read2_strand = "-" if alignment.mate_is_reverse else "+"
                    try:
                        read2_mappings = vci_file.find_mappings(
                            read2_chr, read2_start, read2_end
                        )
                    except:
                        read2_mappings = None

                if read1_mappings is None and read2_mappings is None:

                    alignment_new.flag |= FLAG_UNMAP
                    alignment_new.flag |= FLAG_MUNMAP

                    logger.debug("\tFail due to no mappings")
                    new_file_unmapped.write(alignment)
                    map_statistics_pair["fail_cannot_map"] += 1

                elif read1_mappings is None and read2_mappings and len(read2_mappings) == 1:

                    alignment_new.flag |= FLAG_UNMAP

                    alignment_new.pos = 0
                    alignment_new.cigarstring = "0M"
                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0

                    logger.debug(f"\tPair Success (1:fail,2:simple): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_fail_2_simple"] += 1

                elif read1_mappings is None and read2_mappings and len(read2_mappings) > 1:

                    alignment_new.flag |= FLAG_UNMAP

                    alignment_new.pos = 0
                    alignment_new.cigarstring = "0M"
                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0

                    logger.debug(f"\tPair Success (1:fail,2:complex): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_fail_2_complex"] += 1

                elif read1_mappings and len(read1_mappings) == 1 and read2_mappings is None:

                    alignment_new.flag |= FLAG_MUNMAP

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = alignment.cigar

                    alignment_new.rnext = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pnext = 0
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:simple,2:fail): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_simple_2_fail"] += 1

                elif read1_mappings and len(read1_mappings) == 1 and read2_mappings and len(read2_mappings) == 1:

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = alignment.cigar

                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:simple,2:simple): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_simple_2_simple"] += 1

                elif read1_mappings and len(read1_mappings) == 1 and read2_mappings and len(read2_mappings) > 1:

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = alignment.cigar

                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:simple,2:complex): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_simple_2_complex"] += 1

                elif read1_mappings and len(read1_mappings) > 1 and read2_mappings is None:

                    alignment_new.flag |= FLAG_MUNMAP

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = convert_cigar(alignment.cigar, read_chr, vci_file, alignment.seq, read1_strand, alignment.pos)

                    alignment_new.rnext = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pnext = 0
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:complex,2:fail): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_complex_2_fail"] += 1

                elif read1_mappings and len(read1_mappings) > 1 and read2_mappings and len(read2_mappings) == 1:

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = convert_cigar(alignment.cigar, read_chr, vci_file, alignment.seq, read1_strand, alignment.pos)

                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:complex,2:simple): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_complex_2_simple"] += 1

                elif read1_mappings and len(read1_mappings) > 1 and read2_mappings and len(read2_mappings) > 1:

                    alignment_new.tid = name_to_id[read1_mappings[0].to_chr]
                    alignment_new.pos = read1_mappings[0].to_start
                    alignment_new.cigar = convert_cigar(alignment.cigar, read_chr, vci_file, alignment.seq, read1_strand, alignment.pos)

                    alignment_new.rnext = name_to_id[read2_mappings[0].to_chr]
                    alignment_new.pnext = read2_mappings[0].to_start
                    alignment_new.tlen = 0    # CHECK

                    logger.debug(f"\tPair Success (1:complex,2:complex): {alignment_new.pos} {alignment_new.cigarstring}")
                    new_file.write(alignment_new)
                    map_statistics_pair["success_1_complex_2_complex"] += 1

                else:
                    raise G2GBAMError("Unknown BAM/SAM conversion/parse situation")

    except StopIteration:
        logger.info("All reads processed")

    logger.info(f"  {total:>10} TOTAL ENTRIES")
    logger.info(f"  {total_unmapped:>10} TOTAL UNMAPPED")
    logger.info(f"  {total_fail_qc:>10} TOTAL FAIL QC")

    if map_statistics["total"] > 0:
        logger.info("")
        logger.info("Mapping Summary Single End")
        logger.info(f"  {map_statistics['total']:>10} TOTAL ENTRIES")
        logger.info("")
        logger.info("  {:>10} TOTAL SUCCESS".format(map_statistics["success_simple"] + map_statistics["success_complex"]))
        logger.info("  {:>10} Simple".format(map_statistics["success_simple"]))
        logger.info("  {:>10} Complex".format(map_statistics["success_complex"]))
        logger.info("")
        logger.info("  {:>10} TOTAL FAILURES".format(map_statistics["fail_cannot_map"]))
        logger.info("  {:>10} Cannot Map ".format(map_statistics["fail_cannot_map"]))

    if map_statistics_pair["total"] > 0:
        total_success = 0
        for k, v in map_statistics_pair.items():
            if k.startswith("success"):
                total_success += v

        logger.info("")
        logger.info("Mapping Summary Paired End")
        logger.info("  {:>10} TOTAL ENTRIES".format(map_statistics_pair["total"]))
        logger.info("")
        logger.info("  {:>10} TOTAL SUCCESS".format(total_success))
        logger.info("  {:>10} Read 1 Failed, Read 2 Simple".format(map_statistics_pair["success_1_fail_2_simple"]))
        logger.info("  {:>10} Read 1 Failed, Read 2 Complex".format(map_statistics_pair["success_1_fail_2_complex"]))
        logger.info("  {:>10} Read 1 Simple, Read 2 Failed".format(map_statistics_pair["success_1_simple_2_fail"]))
        logger.info("  {:>10} Read 1 Simple, Read 2 Simple".format(map_statistics_pair["success_1_simple_2_simple"]))
        logger.info("  {:>10} Read 1 Simple, Read 2 Complex".format(map_statistics_pair["success_1_simple_2_complex"]))
        logger.info("  {:>10} Read 1 Complex, Read 2 Failed".format(map_statistics_pair["success_1_complex_2_fail"]))
        logger.info("  {:>10} Read 1 Complex, Read 2 Simple".format(map_statistics_pair["success_1_complex_2_simple"]))
        logger.info("  {:>10} Read 1 Complex, Read 2 Complex".format(map_statistics_pair["success_1_complex_2_complex"]))
        logger.info("")
        logger.info("  {:>10} TOTAL FAILURES".format(map_statistics_pair["fail_cannot_map"]))
        logger.info("  {:>10} Cannot Map".format(map_statistics_pair["fail_cannot_map"]))
        logger.info("")

    logger.info("BAM File Converted")


#
# Functions dealing with CIGAR strings
#
#
#    BAM	OP	    Description
#    0	    M	    alignment match
#    1	    I	    insertion to reference
#    2	    D	    deletion from reference. region deleted from reference genome
#    3	    N	    skipped region from the reference
#    4	    S	    soft clipping (clipped sequence present in SEQ)
#    5	    H	    hard clipping (clipped sequences NOT present in SEQ)
#    6	    P	    padding (silent deletion from padded reference)
#    7	    =	    sequence match
#    8	    X	    sequence mismatch
#

def cigarlist_to_cigarstring(cigar_list: list[tuple(int, int)]) -> str:
    """
    Convert a list of tuples into a cigar string.

    Example::

             [ (0, 10), (1, 1), (0, 75), (2, 2), (0, 20) ]
        =>       10M      1I      75M      2D      20M
        =>                10M1I75M2D20M


    Args:
        cigar_list: A list of tuples where each tuple is a (code, length)

    Returns:
        The cigar string.

    Raises:
        G2GCigarFormatError: On invalid cigar string.
    """
    cigar = ""
    if isinstance(cigar_list, Cigar):
        try:
            for i in cigar_list:
                cigar += str(i.length) + i.code
        except KeyError:
            raise G2GCigarFormatError(f"Invalid cigar code: {i}")
    else:
        try:
            for i in cigar_list:
                cigar += str(i[1]) + CIGAR_N2C[i[0]]
        except KeyError:
            raise G2GCigarFormatError(f"Invalid cigar code: {i}")

    return cigar

def cigar_to_string(cigar):
    """
    Convert a list of tuples into a cigar string.

    Example::

             [ (0, 10), (1, 1), (0, 75), (2, 2), (0, 20) ]
        =>       10M      1I      75M      2D      20M
        =>                10M1I75M2D20M


    :param cigar_list: a list of tuples (code, length)
    :type cigar_list: list
    :return: the cigar string
    :rtype: string
    :raises: :class:`.G2GCigarFormatError` on invalid cigar string
    """
    cigar = ""
    try:
        for i in cigar:
            cigar += str(i.length) + i.code
    except KeyError:
        raise G2GCigarFormatError("Invalid cigar code: " + str(i))

    return cigar

def _cigar_to_list(cigar_string):
    """
    Convert a list of tuples into a cigar string

    Example::

                         10M1I75M2D20M
        =>       10M      1I      75M      2D      20M
        =>   [ (0, 10), (1, 1), (0, 75), (2, 2), (0, 20) ]


    :param cigar_string: a cigar string
    :return: a list of tuples (code, length)
    :rtype: list
    :raises: :class:`.G2GCigarFormatError` on invalid cigar string
    """
    matches = REGEX_CIGAR.findall(cigar_string)
    possible_length = len(REGEX_CIGAR_LENGTH.findall(cigar_string))

    if len(matches) != possible_length:
        raise G2GCigarFormatError("Invalid cigar string: {cigar_string}")

    lst = []
    try:
        for m in matches:
            lst.append(1) # (CIGAR_CODES_REV[m[1]], int(m[0])))
    except KeyError:
        raise G2GCigarFormatError("Invalid cigar string: {0} : {1} ".format(cigar_string, str(m)))

    return lst


def _cigar_convert(cigar, chromosome, vci_file, strand="+", position=0):
    """
    PHASE 1

    Convert each CIGAR element to new mappings and construct an array on NEW cigar elements

    For example, depending on the Intervals in the CHAIN file, let's say we have the following
    CIGAR string: 35M49N65M

    This could get converted into
    35M ==> 4M150D31M
    49N ==> -1N (remember, surrounding M's are used to find the length of N which is done on next pass)
    65M ==> 65M

    First pass yields: 35M49N65M => 4M150D31M-1N65M

    :param cigar:
    :param chromosome:
    :param vci_file:
    :param strand:
    :param position:
    :return:
    """
    cigar_new = []
    current_pos = position
    cigar_no = 0

    for c in cigar:
        cigar_no += 1

        # logger.debug("Element #{0}, '{1}{2}' specified, location: {3}".format(cigar_no, c[1], CIGAR_N2C[c[0]], current_pos))

        increment = c[1]

        if c[0] == CIGAR_m:
            new_mappings = vci_file.find_mappings(chromosome, current_pos, current_pos + c[1])

            if not new_mappings:
                # logger.debug("Mappings: None")
                cigar_new.append(Cigar(CIGAR_S, c[1], 0, 0))
            elif len(new_mappings) == 1:
                # logger.debug("Mappings: Easy: {0}".format(new_mappings[0]))
                cigar_new.append(Cigar(CIGAR_M, new_mappings[0].to_end - new_mappings[0].to_start, new_mappings[0].to_start, new_mappings[0].to_end))
            else:
                # multiple maps, not so easy
                last = None

                for m in new_mappings:
                    # logger.debug("Mappings: Multiple: {0}".format(m))
                    if not last:
                        last = m

                        if current_pos < m.from_start:
                            # special case of first match not in interval, handle accordingly
                            # logger.debug("Adding 'S', because {0} < {1}".format(current_pos, m.from_start))
                            cigar_new.append(Cigar(CIGAR_S, m.from_start - current_pos, 0, 0))
                    else:
                        if m.from_start != last.from_end:
                            # logger.debug("Adding 'M' and 'I', because {0} != {1}".format(m.from_start, last.from_end))

                            cigar_new.append(Cigar(CIGAR_M, last.to_end - last.to_start, last.to_start, last.to_end))
                            cigar_new.append(Cigar(CIGAR_I, m.from_start - last.from_end, last.to_start, last.to_end))

                        elif m.to_start != last.to_end:
                            # logger.debug("Adding 'M' and 'D', because {0} != {1}".format(m.to_start, last.to_end))

                            cigar_new.append(Cigar(CIGAR_M, last.to_end - last.to_start, last.to_start, last.to_end))
                            cigar_new.append(Cigar(CIGAR_D, m.to_start - last.to_end, 0, 0))

                        last = m

                # logger.debug("Adding 'M'")
                cigar_new.append(Cigar(CIGAR_M, last.to_end - last.to_start, last.to_start, last.to_end))

        elif c[0] == CIGAR_i:

            # logger.debug("Adding 'I' and 'D'")

            cigar_new.append(Cigar(CIGAR_I, c[1], 0, 0))
            cigar_new.append(Cigar(CIGAR_D, -1, 0, 0))

            increment = 0

        elif c[0] == CIGAR_d:

            # logger.debug("Adding 'D'")

            cigar_new.append(Cigar(CIGAR_D, -1, 0, 0))

        elif c[0] == CIGAR_n:

            # logger.debug("Adding 'N'")
            cigar_new.append(Cigar(CIGAR_N, -1, 0, 0))

        elif c[0] in [CIGAR_s, CIGAR_h]:

            # logger.debug("Adding '{0}'".format(CIGAR_N2C[c[0]]))
            cigar_new.append(Cigar(CIGAR_N2C[c[0]], c[1], 0, 0))

        else:
            # other
            # logger.debug("OTHER CODE '{0}' found, looking at {1} at {2}".format(CIGAR_N2C[c[0]], c, current_pos))
            raise G2GCigarFormatError("ERROR: Not handling the values in this cigar string: {0}".format(cigar))

        #current_pos += c[1]
        current_pos += increment

        # logger.debug("Current CIGAR: {0}".format(cigar_new))

    return cigar_new


def _cigar_combine_consecutive(cigar):
    """
    Combine consecutive features in a cigar string.

    For example, 2 N's become 1

    :param cigar:
    :return:
    """
    done = False

    while not done:
        done = True

        for i in range(0, len(cigar)-1):

            # logger.debug("{0}={1}".format(i, cigar[i]))
            # logger.debug("{0}={1}".format(i+1, cigar[i+1]))
            if cigar[i].code == cigar[i+1].code:
                done = False
                break

        if not done:
            cigar_temp = []
            cigar_temp.extend(cigar[:i])

            cm1 = cigar[i]
            cm2 = cigar[i+1]
            cm_new = Cigar(cm1.code, cm1.length + cm2.length, cm1.start, cm2.end)
            cigar_temp.append(cm_new)

            # logger.debug("Found consecutive elements {0} and {1}, combined into {2}".format(cm1, cm2, cm_new))

            cigar_temp.extend(cigar[i+2:])

            cigar = cigar_temp

    return cigar


def _cigar_fix_pre_and_post_M(cigar):
    """

    :param cigar:
    :return:
    """
    # pre M to S fix
    for i in range(0, len(cigar)):
        if cigar[i].code == CIGAR_M:
            break

    if i != 0:
        first_m = i

        length = 0
        for i in range(0, first_m):
            if cigar[i].code in [CIGAR_I, CIGAR_S, CIGAR_H]:
                length += cigar[i].length

        temp_cigar = [Cigar(CIGAR_S, length, 0, 0)]
        temp_cigar.extend(cigar[i+1:])

        cigar = temp_cigar

    # post M to S fix
    for i in reversed(range(0, len(cigar))):
        if cigar[i].code == CIGAR_M:
            break

    if i > 0 and i != (len(cigar) - 1):
        last_m = i

        length = 0
        for i in range(last_m+1, len(cigar)):
            if cigar[i].code in [CIGAR_M, CIGAR_I, CIGAR_S]:
                length += cigar[i].length

        temp_cigar = []
        temp_cigar.extend(cigar[:i-1])
        temp_cigar.append(Cigar(CIGAR_S, length, 0, 0))
        cigar = temp_cigar

    return cigar


def _cigar_remove_softs_between_m(cigar):
    """
    Remove soft if surrounded by Ms

    :param cigar:
    :return:
    """
    done = False

    while not done:
        done = True

        for i in range(1, len(cigar)-1):
            if cigar[i].code == CIGAR_S:
                done = False
                break

        if done:
            break

        before = None
        after = None
        for x in reversed(range(i)):
            if cigar[x].code == CIGAR_M:
                before = cigar[x]
                break

        for x in range(i+1, len(cigar)):
            if cigar[x].code == CIGAR_M:
                after = cigar[x]
                break

        if before and after:
            # logger.debug("Found 'S' between 'M' so removing 'S'")
            cigar_temp = []
            cigar_temp.extend(cigar[:i])
            cigar_temp.extend(cigar[i+1:])
            cigar = cigar_temp
            # logger.debug(cigar)
        else:
            done = True

    return cigar


def _cigar_fix_lengths(cigar, sequence):
    """

    :return:
    """
    # Assign length to -1's
    #
    # Since N's aren't mapped we look at the surrounding M's to find the length of the N's
    #
    # Example: 35M49N65M ==> 4M150D31M-1N65M, the -1 will be corrected by finding the last position of the previous
    # M and first position of the next M
    #
    # there are a few special cases that are handled
    # since there were multiple mappings, we will need to figure out the location on the N's
    done = False

    while not done:
        done = True

        # find first element without a length
        i = 0
        for cm in cigar:
            if cm.length == -1:
                break
            i += 1

        if i == len(cigar):
            done = True
            break

        # logger.debug("Found '{0}' at {1}: {2}".format(cm.code, i, cm))

        before = None
        after = None

        # Simple case is surrounded by mapping positions, but might not be the case
        for x in reversed(range(i)):
            if cigar[x].code == CIGAR_M:
                before = cigar[x]
                break


        for x in range(i+1, len(cigar)):
            if cigar[x].code == CIGAR_M:
                after = cigar[x]
                break

        # special case of 89M2000N11M
        # what happens when thi sis converted to 89M-1N11S (no M at end)
        # we should have 89M11S

        # logger.debug("Before: {0}".format(before))
        # logger.debug("After: {0}".format(after))

        # check if all cigar elements from here to end do not have a length
        a = i
        while a < len(cigar) - 1:
            if cigar[a].length != -1:
                break
            a += 1

        # if a == len(cigar_mapping) -1 than all the rest have no length
        # logger.debug("a={0}, len(cigar_mapping) - 1={1}".format(a, len(cigar) - 1))
        if (a == len(cigar) - 1 and cigar[a].start == -1) or not after or not before:
            # take the rest as a clip
            # logger.debug("Found a clip")
            temp_cigar_mappings = cigar[:i]
            temp_total = 0
            for t in temp_cigar_mappings:
                if t.code in [CIGAR_M, CIGAR_I, CIGAR_S]:
                    temp_total += t.length

            temp_cigar_mappings.append(Cigar(CIGAR_S, len(sequence) - temp_total, -1, -1))
            cigar = temp_cigar_mappings
            done = True
        else:
            c = cigar[i]
            new_c = Cigar(c.code, after.start - before.end, before.end, after.start)
            # logger.debug("Replacing, old = {0}, new = {1}".format(c, new_c))
            cigar[i] = new_c

            done = False

    # logger.debug("Removing 0 length elements, if any")
    new_cigar = []
    for cm in cigar:
        if cm[1] == 0:
            # logger.debug("Removing {}".format(cm))
            continue
        new_cigar.append(cm)

    return new_cigar


def convert_cigar(cigar, chromosome, vci_file, sequence, strand='+', position=0, debug_level=0):
    """
    Generate the cigar string of an old alignment.

    P1: Map M with bx; Inherit S and H; Inherit I but put -1D right behind it; Put -1D or -1N when it’s there.
        P1a: Convert xM, if it has zero length after mapping, to xS

    P2: Remove S (including new S originate from unmapped M) if it is surrounded by any pair of consecutive Ms that
    survived P2

    P3: Adjust the size of D or N that are inbetween Ms. Remove it if they have zero length.

    P4: Combine duplicated entries (I guess mostly M or S)

    P5: Put yS for the unmapped regions before the first M and/or after the last M (I believe adding S, H, I’s in
    those regions should get you y). In this phase remove the remaining -1D or -1N in those regions first.

    :param old_alignment: the old alignment
    :type old_alignment: :class:`pysam.AlignedRead`
    :param chromosome:  the chromosome
    :type chromosome: string
    :param vci_filevci_file: the vci_file file
    :type chain: the vci_file file
    :return: a new cigar string based upon the mappings
    :raises: :class:`.G2GCigarFormatError` on invalid cigar string
    """

    old_cigar = cigarlist_to_cigarstring(cigar)
    logger = g2g.get_logger(debug_level)
    logger.debug("CIGAR CONVERSION : {0}".format(old_cigar))

    #
    # PHASE 1: Convert each CIGAR element to new mappings and construct an array on NEW cigar elements
    #

    logger.debug("CIGAR CONVERSION : PHASE 1 : Converting cigar elements")
    new_cigar = _cigar_convert(cigar, chromosome, vci_file, strand, position)
    logger.debug("AFTER PHASE 1 : {0} ".format(new_cigar))

    if len(new_cigar) == 1:

        logger.debug("CIGAR CONVERSION : Skipping to end since only 1 element")

    else:

        #
        # PHASE 2: Remove S if surrounded by M
        #

        logger.debug("CIGAR CONVERSION : PHASE 2 : Remove S if surrounded by M")
        new_cigar = _cigar_remove_softs_between_m(new_cigar)
        logger.debug("AFTER PHASE 2 : {0} ".format(new_cigar))

        #
        # PHASE 3: Fix element lengths
        #

        logger.debug("CIGAR CONVERSION : PHASE 3 : Fix element lengths")
        new_cigar = _cigar_fix_lengths(new_cigar, sequence)
        logger.debug("AFTER PHASE 3 : {0} ".format(new_cigar))

        #
        # PHASE 4: Combine consecutive matching elements
        #

        logger.debug("CIGAR CONVERSION : PHASE 4 : Combining elements")
        new_cigar = _cigar_combine_consecutive(new_cigar)
        logger.debug("AFTER PHASE 4 : {0} ".format(new_cigar))

        #
        # PHASE 5: Combine consecutive matching elements
        #

        logger.debug("CIGAR CONVERSION : PHASE 5 : Fix pre and post Ms")
        new_cigar = _cigar_fix_pre_and_post_M(new_cigar)
        logger.debug("AFTER PHASE 5 : {0} ".format(new_cigar))

    #
    # Final pass through CIGAR string
    #
    # test cigar string length
    #
    # SEQ: segment SEQuence. This field can be a '*' when the sequence is not stored. If not a '*',
    # the length of the sequence must equal the sum of lengths of M/I/S/=/X operations in CIGAR.
    # An '=' denotes the base is identical to the reference base. No assumptions can be made on the
    # letter cases.
    #

    logger.debug("CIGAR CONVERSION : PHASE 6 : Testing length and conversion")

    cigar_seq_length = 0

    # simplify the cigar, throw away the other stuff we used
    simple_cigar = []
    for c in new_cigar:
        simple_cigar.append((CIGAR_C2N[c.code], c.length))
        if c.code in [CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_E, CIGAR_X]:
            cigar_seq_length += c.length

    try:
        if cigar_seq_length != len(sequence):
            logger.debug(
                f"CIGAR SEQ LENGTH={cigar_seq_length} != SEQ_LEN={len(sequence)}"
            )
            # not equal according to chain file format, add the clipping length
            simple_cigar.append((CIGAR_s, len(sequence) - cigar_seq_length))
    except TypeError:
        # avoids: TypeError: object of type 'NoneType' has no len()
        pass

    if old_cigar != cigar_to_string(simple_cigar):
        logger.debug("old cigar != new cigar")
    else:
        logger.debug("old cigar == new cigar")

    logger.debug(f"CIGAR CONVERSION : {old_cigar} ==> {cigar_to_string(simple_cigar)}")

    logger.debug(simple_cigar)
    return simple_cigar


if __name__ == '__main__':
    logger = g2g.get_logger(10)
    cigarstring = '5I3D4M9D3S104M7D2I'
    cigarlist = _cigar_to_list(cigarstring)
    logger.debug(cigarstring)
    logger.debug(cigarlist)
    cigar_new = _cigar_remove_softs_between_m(cigarlist)
    #cigar_new = _cigar_fix_pre_and_post_M(cigarlist)
    print(cigar_to_string(cigar_new))
    print(cigar_new)

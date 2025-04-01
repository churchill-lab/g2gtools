#
# Collection of functions related to VCF files
#
# 1 based
#
#
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
# chr1	2	rs123456	A	G	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	0/0:99:10    # Homozygous reference
# chr1	3	rs123457	A	G,T	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	1/2:99:10    # Heterozygous with two different alt alleles
# chr1	5	rs123458	C	T	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	0/1:99:10    # Heterozygous (ref and first alt)
# chr1	6	rs123459	C	G	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	1/1:99:10    # Homozygous alternate
# chr1	8	rs123460	T	A,C,G	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	0/3:99:10    # Heterozygous (ref and third alt)
# chr1	9	rs123461	T	G	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	./.:99:10    # Missing genotype
# chr1	10	rs123462	G	A	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	1/1:99:10    # Homozygous alternate
# chr1	11	rs123463	G	C	100	PASS	NS=3;DP=30;AF=0.5	GT:GQ:DP	0/1:99:10    # Heterozygous
#
# Position 2: Homozygous reference (0/0) - No change to the reference sequence
# Position 3: Heterozygous with two different alt alleles (1/2) - G in one haplotype, T in the other
# Position 5: Heterozygous with reference and first alt (0/1) - C in one haplotype, T in the other
# Position 6: Homozygous alternate (1/1) - G in both haplotypes
# Position 8: Heterozygous with reference and third alt (0/3) - T in one haplotype, G in the other
# Position 9: Missing genotype (./.) - Typically kept as reference
# Position 10: Homozygous alternate (1/1) - A in both haplotypes
# Position 11: Heterozygous (0/1) - G in one haplotype, C in the other
#
# >EXAMPLE1
# AAACCCTTTGGG
#
# The resulting haplotypes after applying these variants would be:
#
# >EXAMPLE1_haplotype1
# AAGCCTTAGAG
# >EXAMPLE1_haplotype2
# AATCTGTATAC

# standard library imports
from collections import namedtuple
from typing import TextIO
import re

# 3rd party library imports
import pysam

# local library imports
from g2gtools.exceptions import G2GVCFError
from g2gtools.exceptions import G2GValueError
from g2gtools import g2g_utils

VCF_FIELDS = [
    'chrom',
    'pos',
    'id',
    'ref',
    'alt',
    'qual',
    'filter',
    'info',
    'format',
    'samples',
]
VCFRecord = namedtuple('VCFRecord', VCF_FIELDS)

GT_DATA_FIELDS = [
    'ref',
    'left',
    'right',
    'gt',
    'fi',
    'phase',
    'gt_left',
    'gt_right',
    'is_snp',
]
GTData = namedtuple('GTData', GT_DATA_FIELDS)

GENOTYPE_UNPHASED = '/'
GENOTYPE_PHASED = '|'

REGEX_ALT = re.compile(r'(^[A|C|G|T]+)')
REGEX_ALT = re.compile(r'^([ACGT]+)')


class VCFFile:
    """
    Simple VCF object for parsing VCF files.

    This class provides functionality to read and iterate through VCF files,
    supporting transparent gzip decompression. It parses the header to extract
    sample information and provides methods to access genotype data.

    Attributes:
        file_name (str): Path to the VCF file.
        samples (dict[str, int] | None): Dictionary mapping sample names to their column indices.
        current_line (str | None): The current line being processed.
        current_record (VCFRecord | None): The current parsed VCF record.
        reader (TextIO): File reader object for the VCF file.
    """
    file_name: str
    samples: dict[str, int] | None
    current_line: str | None
    current_record: VCFRecord | None
    reader: TextIO

    def __init__(self, file_name: str) -> None:
        """
        Initialize a new VCF file reader.

        Args:
            file_name: The name of the VCF file.

        Raises:
            G2GVCFError: If the VCF header cannot be parsed correctly.
        """
        self.file_name = file_name
        self.samples = None
        self.current_line = None
        self.current_record = None
        self.reader = g2g_utils.open_resource(file_name)
        self.parse_header()

    def parse_header(self) -> None:
        """
        Parse the VCF file header.

        Reads through the metadata lines (starting with ##) and extracts
        sample names from the header line (starting with #).

        Raises:
            G2GVCFError: When the VCF file isn't formatted correctly.
        """
        self.current_line = next(self.reader)

        # Skip metadata lines
        while self.current_line.startswith('##'):
            self.current_line = next(self.reader)

        # Parse column header line
        if self.current_line.startswith('#'):
            elem = self.current_line.strip().split('\t')
            samples = elem[9:]  # Sample names start at column 10

            # Create mapping from sample names to their indices
            self.samples = dict(zip(samples, range(len(samples))))
        else:
            raise G2GVCFError('Improperly formatted VCF file')

    def parse_gt(self, sample: str) -> GTData:
        """
        Parse the GT field from the VCF record.

        Args:
            sample: The sample identifier.

        Returns:
            Parsed GT data into a GTData object.

        Raises:
            G2GVCFError: When the sample doesn't exist or is None.
        """
        if sample is None:
            raise G2GVCFError('Sample must contain a value')

        sample_index = self.get_sample_index(sample)
        return parse_gt(self.current_record, sample_index)

    def __iter__(self) -> 'VCFFile':
        """
        Make the VCF object iterable.

        Returns:
            The VCF object itself as an iterator.
        """
        return self

    def __next__(self) -> VCFRecord:
        """
        Get the next VCF record from the file.

        Returns:
            A VCFRecord object representing the next variant in the file.

        Raises:
            StopIteration: When the end of the file is reached.
        """
        self.current_line = next(self.reader)

        # Skip any remaining header lines
        while self.current_line.startswith('#'):
            self.current_line = next(self.reader)

        # Parse the VCF line into a record
        self.current_record = parse_vcf_line(self.current_line)
        return self.current_record

    # For Python 2 compatibility - in Python 3.10+ this isn't needed
    next = __next__

    def get_sample_index(self, sample: str) -> int:
        """
        Get the sample index.

        Args:
            sample: The sample string.

        Returns:
            The index of the sample in the VCF file.

        Raises:
            G2GVCFError: When the sample doesn't exist.
        """
        if sample in self.samples:
            return self.samples[sample]

        raise G2GVCFError(f"Unknown sample: '{sample}'")

    def close(self) -> None:
        """
        Close the VCF file reader.

        This method should be called when done with the VCF file to release resources.
        """
        if hasattr(self.reader, 'close'):
            self.reader.close()


def parse_vcf_line(line: str | None) -> VCFRecord | None:
    """
    Parse a line in the VCF file.

    Args:
        line: A line (record) from the VCF file.

    Returns:
        A VCFRecord object or None if line cannot be parsed.
    """
    elem = None

    if isinstance(line, str):
        if line.startswith('#'):
            return None

        elem = line.strip().split('\t')
    elif isinstance(line, list):
        elem = line

    try:
        quality = int(elem[5])
    except ValueError:
        try:
            quality = float(elem[5])
        except ValueError:
            quality = None

    filter_field = None

    if elem[6] != '.':
        filter_field = elem[6].split(';')

    info = elem[7]

    try:
        fmt = elem[8]
    except IndexError:
        fmt = None
    else:
        if fmt == '.':
            fmt = None

    return VCFRecord(
        elem[0],
        int(elem[1]),
        None if elem[2] == '.' else elem[2],
        elem[3],
        elem[4].split(','),
        quality,
        filter_field,
        info,
        fmt,
        elem[9:],
    )


def parse_gt(vcf_record: VCFRecord, sample_name: str) -> GTData:
    """
    Parse the GT field from the VCF record.

    Args:
        vcf_record: The current VCF record containing the GT data.
        sample_name:  The name of the sample.

    Returns:
        Parsed GT data into a GTData object.

    Raises:
        G2GVCFError: When an improper VCF field.
        ValueError: When an improper field in the VCF file is bad.
    """
    if sample_name is None:
        raise G2GVCFError('Sample name must contain a value')

    sample_data = vcf_record.samples[sample_name]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find('<') == -1 and sample_data != '.':
        gt_index = vcf_record.format.split(':').index('GT')
        fi_index = vcf_record.format.split(':').index('FI')

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_PHASED)))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_UNPHASED)))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError(f'Unknown phase in GT, {gt}')

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt[genotypes[0] - 1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt[genotypes[1] - 1]

                gt_left = genotypes[0]
                gt_right = genotypes[1]

                # check for to see if ALT is <CN*> or something not ACGT
                if not REGEX_ALT.match(str(gt_left)):
                    left = None
                    gt_left = None

                if not REGEX_ALT.match(str(gt_right)):
                    right = None
                    gt_right = None

        except ValueError:
            pass
        except IndexError:
            pass
        try:
            fi = sample_data.split(':')[fi_index]
        except ValueError:
            pass
        except IndexError:
            pass

    is_snp = (
        len(vcf_record.ref)
        == 1
        == (len(left) if left else 0)
        == (len(right) if right else 0)
    )
    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )


def parse_gt_tuple(
    vcf_record: pysam.VariantRecord, sample_name: str
) -> GTData:
    """
    Parse the GT field within the VCF line.

    Args:
        vcf_record: The VCF record object.
        sample_name: The sample identifier.

    Returns:
        The GT data parsed into a GTData object.

    Raises:
        G2GValueError: When we cannot parse the GT field.
    """
    sample_data = vcf_record.samples[sample_name]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    """
    print(f"{vcf_record.alleles=}")
    print(f"{vcf_record.alleles_variant_types=}")
    print(f"{vcf_record.alts=}")
    print(f"{vcf_record.chrom=}")
    print(f"{vcf_record.contig=}")
    print(f"{vcf_record.pos=}")
    print(f"{vcf_record.qual=}")
    print(f"{vcf_record.ref=}")
    print(f"{vcf_record.rid=}")
    print(f"{vcf_record.rlen=}")
    print(f"{vcf_record.samples=}")
    print(f"{vcf_record.start=}")
    print(f"{vcf_record.stop=}")

    for k, v in sample_data.items():
        print(k, "=", v)
    """

    # check for to see if ALT is <CN*> or something not ACGT
    if (vcf_record.alts is not None) and ','.join(vcf_record.alts).find(
        '<'
    ) == -1:
        try:
            # parse the GT field
            gt = sample_data['GT']
            if any([g is None for g in gt]):
                raise ValueError('bad')

            # assuming no triploids for now
            if gt[0] == 0:
                left = vcf_record.ref
            else:
                left = vcf_record.alts[gt[0] - 1]

            if gt[1] == 0:
                right = vcf_record.ref
            else:
                right = vcf_record.alts[gt[1] - 1]

            gt_left = gt[0]
            gt_right = gt[1]

            fi = vcf_record.format.get('FI')
        except ValueError:
            # LOG.debug(ve)
            pass
        except IndexError:
            # LOG.debug(ie)
            pass

    is_snp = (
        len(vcf_record.ref)
        == 1
        == (len(left) if left else 1)
        == (len(right) if right else 1)
    )
    # gt = GTData(ref='G', left=None, right=None, gt='./.', fi='.', phase=None, gt_left=None, gt_right=None, is_snp=False)
    #  gt=GTData(ref='G', left='G', right='C', gt='0/2', fi='0', phase='/', gt_left=0, gt_right=2, is_snp=True)
    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )


def parse_gt_orig(vcf_record: VCFRecord, sample_index: int) -> GTData:
    """
    Parse the GT field from the VCF record.

    Args:
        vcf_record: The current VCF record containing the GT data.
        sample_index:  The index of the sample.

    Returns:
        Parsed GT data into a GTData object.

    Raises:
        G2GVCFError: When an improper VCF field.
        ValueError: When an improper field in the VCF file is bad.
    """
    if sample_index is None:
        raise G2GVCFError('Sample index must contain a value')

    sample_data = vcf_record.samples[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find('<') == -1 and sample_data != '.':
        gt_index = vcf_record.format.split(':').index('GT')
        fi_index = vcf_record.format.split(':').index('FI')

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_PHASED)))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_UNPHASED)))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError(f'Unknown phase in GT, {gt}')

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt[genotypes[0] - 1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt[genotypes[1] - 1]

                gt_left = genotypes[0]
                gt_right = genotypes[1]

                # check for to see if ALT is <CN*> or something not ACGT
                if not REGEX_ALT.match(gt_left):
                    left = None
                    gt_left = None

                if not REGEX_ALT.match(gt_right):
                    right = None
                    gt_right = None

        except ValueError:
            pass
        except IndexError:
            pass
        try:
            fi = sample_data.split(':')[fi_index]
        except ValueError:
            pass
        except IndexError:
            pass

    is_snp = (
        len(vcf_record.ref)
        == 1
        == (len(left) if left else 0)
        == (len(right) if right else 0)
    )
    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )


def parse_gt_tuple_orig(
    vcf_record: pysam.VCFProxy, sample_index: int
) -> GTData:
    """
    Parse the GT field within the VCF line.

    Args:
        vcf_record: The VCF record object.
        sample_index: The index of the sample.

    Returns:
        The GT data parsed into a GTData object.

    Raises:
        G2GValueError: When we cannot parse the GT field.
    """
    sample_data = vcf_record[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find('<') == -1 and sample_data != '.':
        formats = vcf_record.format.split(':')
        gt_index = formats.index('GT')
        fi_index = formats.index('FI') if 'FI' in formats else None

        try:
            # parse the GT field
            gt = sample_data.split(':')[gt_index]

            # make sure a call can be made
            if gt != '.' and gt != './.' and gt != '.|.':
                if GENOTYPE_PHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_PHASED)))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_UNPHASED)))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise G2GValueError(f'Unknown phase in GT, {gt}')

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt.split(',')[genotypes[0] - 1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt.split(',')[genotypes[1] - 1]

                gt_left = genotypes[0]
                gt_right = genotypes[1]

                # check for to see if ALT is <CN*> or something not ACGT
                # if not REGEX_ALT.match(left) or not REGEX_ALT.match(right):
                #    LOG.error("VFC2VCI CN FOUND")
                #    gt = None
                #    fi = None
                #    left = None
                #    right = None
                #    phase = None
                #    gt_left = None
                #    gt_right = None

        except ValueError:
            # LOG.debug(ve)
            pass
        except IndexError:
            # LOG.debug(ie)
            pass
        try:
            if fi_index:
                fi = sample_data.split(':')[fi_index]
        except ValueError:
            # LOG.debug(ve)
            pass
        except IndexError:
            # LOG.debug(ie)
            pass

    # gt = GTData(ref='G', left=None, right=None, gt='./.', fi='.', phase=None, gt_left=None, gt_right=None, is_snp=False)
    #  gt=GTData(ref='G', left='G', right='C', gt='0/2', fi='0', phase='/', gt_left=0, gt_right=2, is_snp=True)

    is_snp = (
        len(vcf_record.ref)
        == 1
        == (len(left) if left else 1)
        == (len(right) if right else 1)
    )

    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )

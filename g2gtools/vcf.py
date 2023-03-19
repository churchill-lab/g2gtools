#
# Collection of functions related to VCF files
#
# 1 based
#

# standard library imports
from collections import namedtuple
from typing import IO
import re

# 3rd party library imports
# none

# local library imports
from .exceptions import G2GVCFError
from .exceptions import G2GValueError
from . import g2g_utils

VCF_FIELDS = [
    "chrom", "pos", "id", "ref", "alt",
    "qual", "filter", "info", "format", "samples"
]
VCFRecord = namedtuple("VCFRecord", VCF_FIELDS)

GT_DATA_FIELDS = [
    "ref", "left", "right", "gt", "fi",
    "phase", "gt_left", "gt_right", "is_snp"
]
GTData = namedtuple("GTData", GT_DATA_FIELDS)

GENOTYPE_UNPHASED = "/"
GENOTYPE_PHASED = "|"

REGEX_ALT = re.compile(r"(^[A|C|G|T]+)")


class VCFFile(object):
    """
    Simple VCF object for parsing VCF files
    """
    def __init__(self, file_name: str):
        """
        Encapsulate VCF file information.

        Args:
            file_name: The name of the VCF file.
        """
        self.file_name: str = file_name
        self.samples = None
        self.current_line: str | None = None
        self.current_record: VCFRecord | None = None
        self.reader: IO = g2g_utils.open_resource(file_name)
        self.parse_header()

    def parse_header(self):
        """
        Parse the VCF file header.

        Raises:
            G2GVCFError: When the VCF file isn't formatted correctly.
        """
        self.current_line = next(self.reader)

        while self.current_line.startswith("##"):
            self.current_line = next(self.reader)

        if self.current_line.startswith("#"):
            elem = self.current_line.strip().split("\t")
            samples = elem[9:]
            self.samples = dict(zip(samples, (x for x in range(len(samples)))))
        else:
            raise G2GVCFError("Improperly formatted VCF file")

    def parse_gt(self, sample: str):
        """
        Parse the GT field from the VCF record.

        Args:
            sample: The sample identifier.

        Returns:
            Parsed GT data into a GTData object.
        """
        if sample is None:
            raise G2GVCFError("Sample must contain a value")

        sample_index = self.get_sample_index(sample)
        return parse_gt(self.current_record, sample_index)

    def __iter__(self):
        """
        Iterable.
        """
        return self

    def next(self):
        """
        Explicitly call next on the reader.

        Returns:
            A VCFRecord.
        """
        self.current_line = next(self.reader)

        while self.current_line.startswith("#"):
            self.current_line = next(self.reader)

        self.current_record = parse_vcf_line(self.current_line)

        return self.current_record

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
        if line.startswith("#"):
            return None

        elem = line.strip().split("\t")
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

    if elem[6] != ".":
        filter_field = elem[6].split(";")

    info = elem[7]

    try:
        fmt = elem[8]
    except IndexError:
        fmt = None
    else:
        if fmt == ".":
            fmt = None

    return VCFRecord(
        elem[0], int(elem[1]), None if elem[2] == "." else elem[2], elem[3],
        elem[4].split(","), quality, filter_field, info, fmt, elem[9:]
    )


def parse_gt(vcf_record: VCFRecord, sample_index: int) -> GTData:
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
        raise G2GVCFError("Sample index must contain a value")

    sample_data = vcf_record.samples[sample_index]
    gt = None
    fi = None
    left = None
    right = None
    phase = None
    gt_left = None
    gt_right = None

    # check for to see if ALT is <CN*> or something not ACGT
    if vcf_record.alt.find("<") == -1 and sample_data != ".":
        gt_index = vcf_record.format.split(":").index("GT")
        fi_index = vcf_record.format.split(":").index("FI")

        try:
            # parse the GT field
            gt = sample_data.split(":")[gt_index]

            # make sure a call can be made
            if gt != "." and gt != "./." and gt != ".|.":
                if GENOTYPE_PHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_PHASED)))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_UNPHASED)))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise ValueError(f"Unknown phase in GT, {gt}")

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt[genotypes[1]-1]

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
            fi = sample_data.split(":")[fi_index]
        except ValueError:
            pass
        except IndexError:
            pass

    is_snp = len(vcf_record.ref) == 1 == (len(left) if left else 0) == (len(right) if right else 0)
    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )


def parse_gt_tuple(vcf_record: VCFRecord, sample_index: int) -> GTData:
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
    if vcf_record.alt.find("<") == -1 and sample_data != ".":
        formats = vcf_record.format.split(":")
        gt_index = formats.index("GT")
        fi_index = formats.index("FI") if "FI" in formats else None

        try:
            # parse the GT field
            gt = sample_data.split(":")[gt_index]

            # make sure a call can be made
            if gt != "." and gt != "./." and gt != ".|.":
                if GENOTYPE_PHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_PHASED)))
                    phase = GENOTYPE_PHASED
                elif GENOTYPE_UNPHASED in gt:
                    genotypes = list(map(int, gt.split(GENOTYPE_UNPHASED)))
                    phase = GENOTYPE_UNPHASED
                else:
                    raise G2GValueError(f"Unknown phase in GT, {gt}")

                # assuming no triploids for now
                if genotypes[0] == 0:
                    left = vcf_record.ref
                else:
                    left = vcf_record.alt.split(",")[genotypes[0]-1]

                if genotypes[1] == 0:
                    right = vcf_record.ref
                else:
                    right = vcf_record.alt.split(",")[genotypes[1]-1]

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
                fi = sample_data.split(":")[fi_index]
        except ValueError:
            # LOG.debug(ve)
            pass
        except IndexError:
            # LOG.debug(ie)
            pass

    is_snp = len(vcf_record.ref) == 1 == (len(left) if left else 0) == (len(right) if right else 0)
    return GTData(
        vcf_record.ref, left, right, gt, fi, phase, gt_left, gt_right, is_snp
    )

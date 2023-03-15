#
# Collection of g2g utility functions and classes
#

# standard library imports
from itertools import zip_longest
from subprocess import Popen
import bz2
import gzip
import os
import re
import random
import shutil
import string
import sys
import tempfile
import time
import urllib

# 3rd party library imports
from natsort import natsorted as _natsorted
import pysam

# local library imports
from .exceptions import G2GValueError


REGEX_LOCATION = re.compile("(\w*)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
REGEX_LOCATION_CHR = re.compile(
    "(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
BASES = re.compile(r"([ATGCYRSWKMBDHVNatgcyrswkmbdhvn]+)")

TRANS = str.maketrans("ATGCYRSWKMBDHVNatgcyrswkmbdhvn",
                      "TACGRYSWMKVHDBNtacgryswmkvhdbn")

def cmp(x, y):
    """
    cmp(x, y) -> integer

    Return negative if x<y, zero if x==y, positive if x>y.
    """
    return (x > y) - (x < y)

def n(b, encoding="utf-8"):
    return b.decode(encoding)


def s(value):
    if isinstance(value, bytes):
        return n(value)

    if isinstance(value, str):
        return value

    return value


def show_error():
    """
    Show system errors
    """
    et, ev, tb = sys.exc_info()

    print("Error Type: {}".format(et))
    print("Error Value: {}".format(ev))
    print(str(tb))
    while tb:
        co = tb.tb_frame.f_code
        filename = str(co.co_filename)
        line_no = str(tb.tb_lineno)
        print("    {}:{}".format(filename, line_no))
        tb = tb.tb_next


def try_int(s):
    """
    Convert to integer if possible.
    """
    try:
        return int(s)
    except:
        return s


def natsort_key(s):
    """
    Used internally to get a tuple by which s is sorted.
    """
    import re
    return map(try_int, re.findall(r"(\d+|\D+)", s))


def natcmp(a, b):
    """
    Natural string comparison, case sensitive.
    """
    return cmp(natsort_key(a), natsort_key(b))


def natcasecmp(a, b):
    """
    Natural string comparison, ignores case.
    """
    return natcmp(a.lower(), b.lower())


def natsort(seq, cmp=natcmp):
    """
    In-place natural string sort.
    """
    seq.sort(cmp)


def natsorted(seq, cmp=natcmp):
    """
    Returns a copy of seq, sorted by natural string sort.
    """
    return _natsorted(seq)


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def format_time(start, end):
    """
    Format length of time between start and end.

    :param start: the start time
    :param end: the end time
    :return: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def parse_chromosome(chrom_file):
    """
    Parse the chromosome for length
    """
    chromosomes = {}

    try:
        if not os.path.exists(chrom_file):
            raise IOError(chrom_file)

        fd = open(chrom_file, "r")
        for line in fd:
            elem = line.strip().split()
            try:
                chromosome = elem[0]
                chromosome_length = int(elem[1])
                chromosomes[chromosome] = chromosome_length
            except ValueError as e1:
                pass  # most likely header line
        fd.close()
    except IOError as e:
        message = "Error parsing chromosome files: {0}".format(e.message)
        print(message)
        return {}

    return chromosomes


def wrap_sequence(sequence, n=60, fillvalue=""):
    args = [iter(sequence)] * n
    for line in zip_longest(fillvalue=fillvalue, *args):
        yield "".join(line + ("\n",))


def write_sequence(sequence, out, n=60):
    for i in range(0, len(sequence), n):
        out.write(sequence[i:i + n])
        out.write("\n")


def dump_file_contents(file_name):
    with open(file_name, "r") as f:
        shutil.copyfileobj(f, sys.stdout)


def reverse_sequence(sequence):
    """
    Return the reverse of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse string or Sequence object
    """
    if isinstance(sequence, str):
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
        raise ValueError(f"Sequence contains non-DNA character '{val[position]}' at position {position+1:n}\n")

    if isinstance(sequence, str):
        return val
    else:
        raise ValueError(f"Cannot complement object of {type(sequence)}, expecting string or Sequence object")


def reverse_complement_sequence(sequence):
    """
    Return the reverse-complement of sequence

    :param sequence: either a string or Sequence object
    :return: the reverse-complement string or Sequence object
    """
    return reverse_sequence(complement_sequence(sequence))


def concatenate_files(
        list_files: list[str],
        to_file: str,
        delete: bool | None = False,
        mode: str | None = "wb"
) -> None:
    """
    Concatenate the files in the order they are listed.

    list_files: A list of the files to concatenate.
    to_file: Name of the file to concatenate to.
    delete: True to delete the files after they have been concatenated.
    mode: Mode to write the concatenated file.
    """
    with open(to_file, mode) as out:
        for from_file in list_files:
            with open(from_file, "rb") as f:
                shutil.copyfileobj(f, out)
            if delete:
                delete_file(from_file)


def delete_file(filename: str) -> None:
    """
    Delete the file.

    Args:
        filename: Name of the file to parse.
    """
    try:
        os.remove(filename)
    except OSError:
        pass


def delete_index_files(filename):
    """
    Delete the file indices of the filename.
        filename.gzi
        filename.fai
        filename.tbi

    Args:
        filename: Name of the file to parse.
    """
    delete_file(f"{filename}.gzi")
    delete_file(f"{filename}.fai")
    delete_file(f"{filename}.tbi")


def bgzip_decompress(filename):
    # gobble the output
    with open(os.devnull, "w") as fnull:
        p = Popen(["bgzip", "-f", "-d", filename], stdout=fnull, stderr=fnull)
        p.wait()

    delete_index_files(filename)


def bgzip_file(
    original_file: str,
    new_file: str,
    delete_original: bool | None = False,
    force: bool | None = True
) -> None:
    """
    bgzip a file and index it

    Args:
        original_file: The file to compress and index.
        new_file: Name of the new file.
        delete_original: True to delete the original file.
        force: True to force overwrite and index.
    """
    pysam.tabix_compress(original_file, new_file, force)

    if delete_original:
        delete_file(original_file)


def has_index_file(original_file, file_format=None):
    """

    :param original_file:
    :param new_file:
    :param file_format:
    :return:
    """

    if not file_format:
        # try to guess the file format
        if original_file.lower().endswith(".fa") or original_file.lower().endswith(".fasta"):
            file_format = "fa"
        elif original_file.lower().endswith(".vcf"):
            file_format = "vcf"
        elif original_file.lower().endswith(".vci"):
            file_format = "vci"
        else:
            raise G2GValueError("Cannot determine file format")

    if file_format.lower() == "fa":
        ext = "fai"
    elif file_format.lower() == "vcf":
        ext = "tbi"
    elif file_format.lower() == "vci":
        ext = "tbi"
    else:
        raise G2GValueError(f"Unknown file format: {file_format}")

    idx_file = f"{original_file}.{ext}"

    return os.path.exists(idx_file)


def index_file(
        original_file: str,
        file_format: str | None = "vcf",
        overwrite: bool | None = False
) -> None:
    """
    Parse the VCF file and create a VCI file.

    Args
        original_file (str): Name of the file to index.
        file_format (str): Format of the file (fa, vcf, vci).
        overwrite(bool): True to overwrite existing file.
    """
    if overwrite or not has_index_file(original_file, file_format=file_format):
        if file_format.lower() == "fa":
            pysam.FastaFile(original_file)
        elif file_format.lower() == "vcf":
            pysam.tabix_index(original_file, preset="vcf", force=True)
        elif file_format.lower() == "vci":
            pysam.tabix_index(
                original_file, seq_col=0, start_col=1, end_col=1, force=True
            )
        else:
            raise G2GValueError(f"Unknown file format: {file_format}")


def bgzip_and_index_file(
        original_file: str,
        new_file: str,
        delete_original: bool | None = False,
        force: bool | None = True,
        file_format: str | None = "vcf"
) -> None:
    """
    bgzip a file and index it

    Args:
        original_file: The file to compress and index.
        new_file: Name of the new file.
        delete_original: True to delete the original file.
        force: True to force overwrite and index.
        file_format: Format of the file so we can compress and index correctly.
    """
    bgzip_file(original_file, new_file, delete_original, force)
    index_file(new_file, file_format)


def open_resource(resource, mode="rb"):
    """
    Open different types of files and return the handle.

    :param resource: a file located locally or on the internet.  Gzip'd and zip'd files are handled.
    :param mode: standard file open modes
    :return: the resource (file) handle
    """
    if not resource:
        return None

    if not isinstance(resource, str):
        return resource

    resource = s(resource)

    if resource.endswith((".gz", ".Z", ".z")):
        return gzip.open(resource, mode)
    elif resource.endswith((".bz", ".bz2", ".bzip2")):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(("http://", "https://", "ftp://")):
        return urllib.urlopen(resource)
    else:
        return open(resource, mode)


def create_random_string(size=6, chars=string.ascii_uppercase + string.digits):
    return "".join(random.choice(chars) for _ in range(size))


def gen_file_name(
        name: str | None = None,
        prefix: str | None = "",
        output_dir=".",
        extension="log",
        append_time=True
) -> str:
    """
    Generate a file name.

    Args
        name: A base name for the file.
        prefix: A prefix for the file if necessary.
        output_dir: name of the directory
        extension: which strain to process
        append_time: True to place troubling VCF lines in extra file

    Returns
        The absolute path to the file.
    """
    if name is None:
        name = create_random_string(15)

    if extension and extension[-1] == ".":
        extension = extension[-1:]

    if append_time:
        t = time.strftime("%Y-%m-%d.%H_%M_%S")
        file_name = f"{prefix}{name}.{t}.{extension}"
    else:
        if extension:
            file_name = f"{prefix}{name}.{extension}"
        else:
            file_name = f"{prefix}{name}"

    return os.path.abspath(os.path.join(output_dir, file_name))


def get_extension(filename):
    file_first, file_extension = os.path.splitext(filename)
    return file_extension


def prepend_before_extension(filename, text):
    file_first, file_extension = os.path.splitext(filename)

    if file_extension:
        return f"{file_first}.{text}{file_extension}"

    return filename


def get_dir_and_file(filename):
    abspath = os.path.abspath(filename)
    if os.path.isdir(abspath):
        return abspath, None

    return os.path.split(abspath)


def check_file(
        file_name: str,
        mode: str | None = "r"
) -> str:
    """
    Check if file_name exists and accessible for reading or writing.

    Args:
        file_name: The name of the file.
        mode: "r" for reading, "w" for writing.

    Returns:
        The absolute path of the file.
    """
    if mode == "r":
        if file_name and os.path.exists(file_name):
            return os.path.abspath(file_name)

        raise G2GValueError(f"The following file does not exist: {file_name}")
    elif mode == "w":
        file_dir = "."

        if file_name:
            file_name = os.path.abspath(file_name)
            file_dir = os.path.dirname(file_name)

            if not os.access(file_dir, os.W_OK | os.X_OK):
                raise G2GValueError(f"Cannot generate file: {file_name}")

            return file_name

    raise G2GValueError(f"Unspecified mode to open file, '{mode}'")


def create_temp_dir(name, prefix=".g2gtools_", dir=None):
    new_name = f"{prefix}{name}"
    return os.path.abspath(tempfile.mkdtemp(prefix=new_name, dir=dir))


def get_sys_exec_root_or_drive():
    path = sys.executable
    while os.path.split(path)[1]:
        path = os.path.split(path)[0]
    return path


def delete_dir(dir):
    dir = os.path.abspath(dir)
    if os.path.exists(dir):
        root = get_sys_exec_root_or_drive()
        # bad error checking, but at least it's something
        if root == dir:
            raise G2GValueError(f"Will not delete directory: {dir}")
        try:
            shutil.rmtree(dir)
        except Exception as e:
            raise G2GValueError(f"Will not delete directory: {dir}")


def location_to_filestring(location):
    return f"{location.seq_id}-{location.start}-{location.end}"


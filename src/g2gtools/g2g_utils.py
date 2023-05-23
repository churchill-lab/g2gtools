# standard library imports
from itertools import zip_longest
from os import PathLike
from subprocess import Popen
from typing import Any
from typing import IO
from typing import Iterator
import bz2
import gzip
import logging
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
from rich.logging import RichHandler
from natsort import natsorted as _natsorted
import pysam

# local library imports
import g2gtools.g2g
from g2gtools.exceptions import G2GValueError


REGEX_LOCATION = re.compile(
    r'(\w*)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?',
    re.IGNORECASE,
)
REGEX_LOCATION_CHR = re.compile(
    r'(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?',
    re.IGNORECASE,
)
BASES = re.compile(r'([ATGCYRSWKMBDHVNatgcyrswkmbdhvn]+)')

TRANS = str.maketrans(
    'ATGCYRSWKMBDHVNatgcyrswkmbdhvn', 'TACGRYSWMKVHDBNtacgryswmkvhdbn'
)


def get_logger(logger_name: str = 'g2gtools') -> logging.Logger:
    """
    Get the logger.

    Args:
        logger_name: Name of the logger.

    Returns:
        logging.Logger: The logging object.
    """
    return logging.getLogger(logger_name)


def configure_logging(
    logger_name: str = 'g2gtools', level: int = 0
) -> logging.Logger:
    """
    Configure the logger with the specified `level`. Valid `level` values
    are:

    ======  =================================
    level   logging value
    ======  =================================
    0       logging.WARNING is informational
    1       logging.INFO is user debug
    2+      logging.DEBUG is developer debug
    ======  =================================

    Anything greater than 2 is treated as 2.

    if the environment variable G2GTOOLS_APP_DEBUG is set to 1,
    there will be more detailed debugging information.

    Args:
        logger_name: The name of the logger.
        level: The logging level; defaults to 0.

    Returns:
        logging.Logger: The logging object.
    """
    ensimpl_app_debug = nvli(os.environ.get('G2GTOOLS_APP_DEBUG', '0'), -1)

    rich_handler = RichHandler(
        level=logging.NOTSET,
        show_level=False,
        show_time=True,
        show_path=False,
        omit_repeated_times=False,
    )

    if ensimpl_app_debug == 1:
        rich_handler = RichHandler(
            level=logging.NOTSET,
            show_level=True,
            show_time=True,
            show_path=True,
            omit_repeated_times=False,
        )

    # this is configuring the root logger and below
    # set the level to WARNING, as that is the default
    logging.basicConfig(
        level=logging.WARNING,
        format='%(message)s',
        datefmt=f'{logger_name} [%X]',
        handlers=[rich_handler],
    )

    log = logging.getLogger(logger_name)

    # set gbrs's logging level
    if level == 0:
        log.setLevel(logging.WARNING)
    elif level == 1:
        log.setLevel(19)
    elif level > 1:
        log.setLevel(logging.DEBUG)

    return log


def nvli(value, default) -> int:
    """Returns `value` as an int if `value` can be converted, else `default`.

    Args:
        value: The value to evaluate and convert to an it.
        default: The default value.

    Returns:
        Either `value` or `default`.
    """
    ret = default
    if value:
        try:
            ret = int(value)
        except ValueError:
            pass
    return ret


def cmp(x: Any, y: Any) -> int:
    """
    Compare x and y.

    Args:
        x: The first value.
        y: The second value.

    Returns:
        Negative if x<y, zero if x==y, positive if x>y.
    """
    return (x > y) - (x < y)


def n(b, encoding='utf-8'):
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

    print('Error Type: {}'.format(et))
    print('Error Value: {}'.format(ev))
    print(str(tb))
    while tb:
        co = tb.tb_frame.f_code
        filename = str(co.co_filename)
        line_no = str(tb.tb_lineno)
        print('    {}:{}'.format(filename, line_no))
        tb = tb.tb_next


def try_int(s: Any) -> Any:
    """
    Convert to integer if possible.

    Args:
        s: The value to convert to an int.

    Returns:
        The value converted to int if possible, else the same value.
    """
    try:
        return int(s)
    except Exception:
        return s


def natsort_key(s):
    """
    Used internally to get a tuple by which s is sorted.
    """
    return map(try_int, re.findall(r'(\d+|\D+)', s))


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


def format_time(start: float, end: float) -> str:
    """
    Format length of time between start and end.

    Args:
        start: The start time.
        end: The end time.

    Returns:
        A formatted string of hours, minutes, and seconds.
    """
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return '{:0>2}:{:0>2}:{:05.2f}'.format(int(hours), int(minutes), seconds)


def parse_chromosome(chrom_file_name: str) -> dict[str, int]:
    """
    Parse the chromosome file for length.

    Args:
        chrom_file_name: Name of the file to parse.

    Returns:
        A dictionary with chromosome as keys and the lengths as values.
    """
    chromosomes = {}

    try:
        if not os.path.exists(chrom_file_name):
            raise IOError(chrom_file_name)

        fd = open(chrom_file_name, 'r')
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
        message = f'Error parsing chromosome files: {e.message}'
        print(message)
        return {}

    return chromosomes


def wrap_sequence(
    sequence: str, wrap_length: int = 60, fill_value: str | None = ''
) -> Iterator[str]:
    """
    Wrap a sequence for better output.

    Args:
        sequence: A sequence string.
        wrap_length: How many characters before wrapping.
        fill_value: Fill in value.

    Yields:
        Iterator of str.
    """
    args = [iter(sequence)] * wrap_length
    for line in zip_longest(fillvalue=fill_value, *args):
        yield ''.join(line + ('\n',))


def write_sequence(sequence: str, out: IO, wrap_length: int = 60):
    """
    Write the formatted sequence to output.

    Args:
        sequence: A sequence string.
        out: Output to write to.
        wrap_length: How many characters before wrapping.
    """
    for i in range(0, len(sequence), wrap_length):
        out.write(sequence[i : i + wrap_length])
        out.write('\n')


def dump_file_contents(file_name: str) -> None:
    """
    Display the contents of the file to stdout.

    Args:
        file_name: The name of the file.
    """
    with open(file_name, 'r') as f:
        shutil.copyfileobj(f, sys.stdout)


def reverse_sequence(sequence: str) -> str:
    """
    Return the reverse of sequence.

    Args:
        sequence: A sequence string.

    Returns:
        The reverse of the sequence string.
    """
    if isinstance(sequence, str):
        return sequence[::-1]
    else:
        raise ValueError(
            f'Cannot complement object of {type(sequence)}, '
            'expecting string or Sequence object'
        )


def complement_sequence(sequence):
    """
    Return the complement of sequence.

    Args:
        sequence: A sequence string.

    Returns:
        The complement of the sequence string.
    """
    val = str(sequence)

    if len(re.findall(BASES, val)) == 1:  # finds invalid characters if > 1
        val = str(sequence).translate(TRANS)
    else:
        matches = re.findall(BASES, val)
        position = len(matches[0])
        raise ValueError(
            'Sequence contains non-DNA character '
            f"'{val[position]}' at position {position+1:n}"
        )

    if isinstance(sequence, str):
        return val
    else:
        raise ValueError(
            f'Cannot complement object of {type(sequence)}, '
            'expecting string or Sequence object'
        )


def reverse_complement_sequence(sequence):
    """
    Return the reverse complement of sequence string.

    Args:
        sequence: A sequence string.

    Returns:
        The reverse complement of the sequence string.
    """
    return reverse_sequence(complement_sequence(sequence))


def concatenate_files(
    list_files: list[str],
    file_name_out: str,
    delete: bool | None = False,
    mode: str | None = 'wb',
) -> None:
    """
    Concatenate the files in the order they are listed.

    list_files: A list of the files to concatenate.
    to_file: Name of the file to concatenate to.
    delete: True to delete the files after they have been concatenated.
    mode: Mode to write the concatenated file.
    """
    with open(file_name_out, mode) as out:
        for from_file in list_files:
            with open(from_file, 'rb') as f:
                shutil.copyfileobj(f, out)
            if delete:
                delete_file(from_file)


def delete_file(file_name: str) -> None:
    """
    Delete the file.

    Args:
        file_name: Name of the file to parse.
    """
    try:
        os.remove(file_name)
    except OSError:
        pass


def delete_index_files(file_name: str) -> None:
    """
    Delete the file indices of the filename.
        filename.gzi
        filename.fai
        filename.tbi

    Args:
        file_name: Name of the file to delete the index of.
    """
    delete_file(f'{file_name}.gzi')
    delete_file(f'{file_name}.fai')
    delete_file(f'{file_name}.tbi')


def bgzip_decompress(file_name: str) -> None:
    """
    Decompress the file specified by file_name.

    Args:
        file_name: The file to bgzip decompress.
    """
    # gobble the output
    with open(os.devnull, 'w') as null:
        p = Popen(['bgzip', '-f', '-d', file_name], stdout=null, stderr=null)
        p.wait()

    delete_index_files(file_name)


def bgzip_file(
    file_name_in: str,
    file_name_out: str,
    delete_original: bool | None = False,
    force: bool | None = True,
) -> None:
    """
    bgzip a file and index it.

    Args:
        file_name_in: The file to compress and index.
        file_name_out: Name of the new file.
        delete_original: True to delete the original file.
        force: True to force overwrite and index.
    """
    pysam.tabix_compress(file_name_in, file_name_out, force)

    if delete_original:
        delete_file(file_name_in)


def has_index_file_fasta(file_name: str) -> bool:
    """
    Check to see if the Fasta file has an index.

    Args:
        file_name: The name of the Fasta file to check.

    Returns:
        True if the index exists, False otherwise.

    Raises:
        G2GValueError: If the index cannot be determined.
    """
    if file_name.lower().endswith('.fa'):
        fai_file = f'{file_name}.fai'
        return os.path.exists(fai_file)
    elif file_name.lower().endswith('.gz'):
        fai_file = f'{file_name}.fai'
        gzi_file = f'{file_name}.gzi'
        return os.path.exists(fai_file) and os.path.exists(gzi_file)

    raise G2GValueError('Cannot determine file format')


def has_index_file_vcf(file_name: str) -> bool:
    """
    Check to see if the VCF file has an index.

    Args:
        file_name: The name of the VCF file to check.

    Returns:
        True if the index exists, False otherwise.

    Raises:
        G2GValueError: If the index cannot be determined.
    """
    if file_name.lower().endswith('.vcf.gz'):
        tbi_file = f'{file_name}.tbi'
        csi_file = f'{file_name}.csi'
        return os.path.exists(tbi_file) or os.path.exists(csi_file)
    elif file_name.lower().endswith('.vcf'):
        raise G2GValueError('VCF files should be compressed')

    raise G2GValueError('Cannot determine file format')


def has_index_file_vci(file_name: str) -> bool:
    """
    Check to see if the VCI file has an index.

    Args:
        file_name: The name of the VCI file to check.

    Returns:
        True if the index exists, False otherwise.

    Raises:
        G2GValueError: If the index cannot be determined.
    """
    if file_name.lower().endswith('.vci.gz'):
        tbi_file = f'{file_name}.tbi'
        csi_file = f'{file_name}.csi'
        return os.path.exists(tbi_file) or os.path.exists(csi_file)
    elif file_name.lower().endswith('.vci'):
        raise G2GValueError('VCI files should be compressed')

    raise G2GValueError('Cannot determine file format')


def has_index_file(file_name: str) -> bool:
    """
    Check to see if the file has an index.

    Args:
        file_name: The name of the file to check.

    Returns:
        True if the index exists, False otherwise.

    Raises:
        G2GValueError: If the index cannot be determined.
    """
    if file_name.lower().endswith('.fa'):
        return has_index_file_fasta(file_name)
    elif file_name.lower().endswith('.fasta'):
        return has_index_file_fasta(file_name)
    elif file_name.lower().endswith('.vcf'):
        return has_index_file_vcf(file_name)
    elif file_name.lower().endswith('.vcf.gz'):
        return has_index_file_vcf(file_name)
    elif file_name.lower().endswith('.vci'):
        return has_index_file_vci(file_name)
    elif file_name.lower().endswith('.vci.gz'):
        return has_index_file_vci(file_name)
    else:
        raise G2GValueError('Cannot determine file format')


def index_file(
    file_name: str, file_format: str = 'vcf', overwrite: bool = False
) -> None:
    """
    Parse the VCF file and create a VCI file.

    Args
        file_name: Name of the file to index.
        file_format: Format of the file (fa, vcf, vci).
        overwrite: True to overwrite existing file.
    """
    if overwrite or not has_index_file(file_name):
        if file_format.lower() == 'fa':
            pysam.FastaFile(file_name)
        elif file_format.lower() == 'vcf':
            pysam.tabix_index(file_name, preset='vcf', force=True)
        elif file_format.lower() == 'vci':
            pysam.tabix_index(
                file_name, seq_col=0, start_col=1, end_col=1, force=True
            )
        else:
            raise G2GValueError(f'Unknown file format: {file_format}')


def bgzip_and_index_file(
    file_name_in: str,
    file_name_out: str,
    delete_original: bool = False,
    force: bool = True,
    file_format: str = 'vcf',
) -> None:
    """
    bgzip a file and index it

    Args:
        file_name_in: The file to compress and index.
        file_name_out: Name of the new file.
        delete_original: True to delete the original file.
        force: True to force overwrite and index.
        file_format: Format of the file, so we can compress and index correctly.
    """
    bgzip_file(file_name_in, file_name_out, delete_original, force)
    index_file(file_name_out, file_format)


def open_resource(resource, mode: str = 'rb'):
    """
    Open different types of files and return the handle.

    Args:
        resource: A file located locally or on the internet.
            Gzip'd and zip'd files are handled.
        mode: Standard file open modes.

    Returns:
        The resource (file) handle.
    """
    if not resource:
        return None

    if not isinstance(resource, str):
        return resource

    resource = s(resource)

    if resource.endswith(('.gz', '.Z', '.z')):
        return gzip.open(resource, mode)
    elif resource.endswith(('.bz', '.bz2', '.bzip2')):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(('http://', 'https://', 'ftp://')):
        return urllib.urlopen(resource)
    else:
        return open(resource, mode)


def create_random_string(
    size: int = 6, chars: str = string.ascii_uppercase + string.digits
) -> str:
    """
    Create a random string.

    Args:
        size: The length of the string.
        chars: The character set to use.

    Returns:
        A random string of length size.
    """
    return ''.join(random.choice(chars) for _ in range(size))


def gen_file_name(
    name: str | None = None,
    prefix: str = '',
    output_dir: str = '.',
    extension: str = 'log',
    append_time: bool = True,
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

    if extension and extension[-1] == '.':
        extension = extension[-1:]

    if append_time:
        t = time.strftime('%Y-%m-%d.%H_%M_%S')
        file_name = f'{prefix}{name}.{t}.{extension}'
    else:
        if extension:
            file_name = f'{prefix}{name}.{extension}'
        else:
            file_name = f'{prefix}{name}'

    return os.path.abspath(os.path.join(output_dir, file_name))


def get_extension(filename):
    file_first, file_extension = os.path.splitext(filename)
    return file_extension


def prepend_before_extension(filename, text):
    file_first, file_extension = os.path.splitext(filename)

    if file_extension:
        return f'{file_first}.{text}{file_extension}'

    return filename


def get_dir_and_file(filename):
    abspath = os.path.abspath(filename)
    if os.path.isdir(abspath):
        return abspath, None

    return os.path.split(abspath)


def check_file(file_name: str, mode: str | None = 'r') -> str:
    """
    Check if file_name exists and accessible for reading or writing.

    Args:
        file_name: The name of the file.
        mode: "r" for reading, "w" for writing.

    Returns:
        The absolute path of the file.

    Raises:
        G2GValueError: When the file doesn't exist or cannot be generated.
    """
    if mode == 'r':
        if file_name and os.path.exists(file_name):
            return os.path.abspath(file_name)

        raise G2GValueError(f'The following file does not exist: {file_name}')
    elif mode == 'w':
        file_dir = '.'

        if file_name:
            file_name = os.path.abspath(file_name)
            file_dir = os.path.dirname(file_name)

            if not os.access(file_dir, os.W_OK | os.X_OK):
                raise G2GValueError(f'Cannot generate file: {file_name}')

            return file_name

    raise G2GValueError(f"Unspecified mode to open file, '{mode}'")


def create_temp_dir(
    name: str, prefix: str = '.g2gtools_', dir: str | PathLike | None = None
) -> str:
    """
    Create a temporary directory.

    Args:
        name: Name of the directory.
        prefix: Prefix of the name.
        dir: Directory where to create temp directory.

    Returns:
        The absolute path to the temporary directory.
    """
    new_name = f'{prefix}{name}'
    return os.path.abspath(tempfile.mkdtemp(prefix=new_name, dir=dir))


def get_sys_exec_root_or_drive() -> str:
    """
    Get system root directory.

    Returns:
        The path to the root dir.
    """
    path = sys.executable
    while os.path.split(path)[1]:
        path = os.path.split(path)[0]
    return path


def delete_dir(directory: str) -> None:
    """
    Delete the specified directory.

    Args:
        directory: The directory to delete

    Raises:
        G2GValueError: When the directory cannot be deleted.
    """
    directory = os.path.abspath(directory)
    if os.path.exists(directory):
        root = get_sys_exec_root_or_drive()
        # bad error checking, but at least it's something
        if root == directory:
            raise G2GValueError(f'Will not delete directory: {directory}')
        try:
            shutil.rmtree(directory)
        except Exception as e:
            raise G2GValueError(f'Will not delete directory: {directory}')


def location_to_filestring(location: g2gtools.g2g.Region) -> str:
    """
    Convert the specified Region into a string representation.

    Args:
        location: The Region object.

    Returns:
        The region as a string.
    """
    return f'{location.seq_id}-{location.start}-{location.end}'

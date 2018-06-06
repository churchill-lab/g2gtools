# -*- coding: utf-8 -*-

#
# Collection of g2g utility functions and classes
#

from __future__ import print_function
from __future__ import with_statement
from past.builtins import cmp, xrange
from future.utils import bytes_to_native_str as n
from future.moves.itertools import zip_longest

import re
import random
import string
import sys
import tempfile
import traceback

from natsort import natsorted as _natsorted

from . import exceptions
from . import compat

REGEX_LOCATION = re.compile("(\w*)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
REGEX_LOCATION_CHR = re.compile("(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)
BASES = re.compile(r"([ATGCYRSWKMBDHVNatgcyrswkmbdhvn]+)")

if compat.is_py2:
    TRANS = string.maketrans('ATGCYRSWKMBDHVNatgcyrswkmbdhvn',
                             'TACGRYSWMKVHDBNtacgryswmkvhdbn')
else:
    TRANS = str.maketrans('ATGCYRSWKMBDHVNatgcyrswkmbdhvn',
                          'TACGRYSWMKVHDBNtacgryswmkvhdbn')



def s(value):
    if compat.is_py2:
        if isinstance(value, unicode):
            return value.decode('ascii', 'ignore')

        if isinstance(value, str):
            return value

    else:
        if isinstance(value, bytes):
            return n(value)

        if isinstance(value, str):
            return value

    return value

def _show_error():
    """
    show system errors
    """
    et, ev, tb = sys.exc_info()

    print("Error Type: {}".format(et))
    print("Error Value: {}".format(ev))
    print(str(tb))
    while tb :
        co = tb.tb_frame.f_code
        filename = str(co.co_filename)
        line_no = str(tb.tb_lineno)
        print('    {}:{}'.format(filename, line_no))
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


def format_time(start, end):
    """
    Format length of time between start and end.

    :param start: the start time
    :param end: the end time
    :return: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def parse_chromosome(chrom_file):
    '''
    Parse the chromosome for length
    '''
    chromosomes = {}

    try:
        if not os.path.exists(chrom_file):
            raise IOError(chrom_file)

        fd = open(chrom_file, 'r')
        for line in fd:
            elem = line.strip().split()
            try:
                chromosome = elem[0]
                chromosome_length = int(elem[1])
                chromosomes[chromosome] = chromosome_length
            except ValueError as e1:
                pass # most likely header line
        fd.close()
    except IOError as e:
        message = "Error parsing chromosome files: {0}".format(e.message)
        print(message)
        return {}

    return chromosomes


def wrap_sequence(sequence, n=60, fillvalue=''):
    args = [iter(sequence)] * n
    for line in zip_longest(fillvalue=fillvalue, *args):
        yield ''.join(line + ('\n',))


def write_sequence(sequence, out, n=60):
    for i in xrange(0, len(sequence), n):
        out.write(sequence[i:i+n])
        out.write('\n')


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
        raise ValueError("Sequence contains non-DNA character '{0}' at position {1:n}\n".format(val[position], position + 1))

    if isinstance(sequence, str):
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



import bz2
import gzip
import os
import shutil
import time
import urllib
from subprocess import Popen

import pysam

from .exceptions import G2GValueError


def concatenate_files(list_files, to_file, delete=False, mode='wb'):
    with open(to_file, mode) as out:
        for from_file in list_files:
            with open(from_file, 'rb') as f:
                shutil.copyfileobj(f, out)
            if delete:
                delete_file(from_file)


def delete_file(file_name):
    """
    :param file_name:
    :return:
    """
    try:
        os.remove(file_name)
    except:
        pass


def delete_index_files(filename):
    delete_file("{0}.gzi".format(filename))
    delete_file("{0}.fai".format(filename))
    delete_file("{0}.tbi".format(filename))


def bgzip_decompress(filename):
    # gobble the output
    with open(os.devnull, "w") as fnull:
        p = Popen(['bgzip', '-f', '-d', filename], stdout=fnull, stderr=fnull)
        p.wait()

    delete_index_files(filename)


def bgzip_file(original_file, new_file, delete_original=False, force=True):
    """

    :param original_file:
    :param new_file:
    :param force:
    :return:
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
            file_format = 'fa'
        elif original_file.lower().endswith(".vcf"):
            file_format = 'vcf'
        elif original_file.lower().endswith(".vci"):
            file_format = 'vci'
        else:
            raise G2GValueError("Cannot determine file format")

    if file_format.lower() == 'fa':
        ext = 'fai'
    elif file_format.lower() == 'vcf':
        ext = 'tbi'
    elif file_format.lower() == 'vci':
        ext = 'tbi'
    else:
        raise G2GValueError("Unknown file format: {0}".format(file_format))

    idx_file = '{}.{}'.format(original_file, ext)

    return os.path.exists(idx_file)


def index_file(original_file, file_format="vcf", overwrite=False):
    """

    :param original_file:
    :param new_file:
    :param file_format:
    :return:
    """
    if overwrite or not has_index_file(original_file, file_format=file_format):
        if file_format.lower() == 'fa':
            pysam.faidx(original_file)
        elif file_format.lower() == 'vcf':
            pysam.tabix_index(original_file, preset="vcf", force=True)
        elif file_format.lower() == 'vci':
            pysam.tabix_index(original_file, seq_col=0, start_col=1, end_col=1, force=True)
        else:
            raise G2GValueError("Unknown file format: {0}".format(file_format))


def bgzip_and_index_file(original_file, new_file, delete_original=False, force=True, file_format="vcf"):
    bgzip_file(original_file, new_file, delete_original, force)
    index_file(new_file, file_format)


def open_resource(resource, mode='rb'):
    """
    Open different types of files and return the handle.

    :param resource: a file located locally or on the internet.  Gzip'd and zip'd files are handled.
    :param mode: standard file open modes
    :return: the resource (file) handle
    """
    if not resource:
        return None

    if compat.is_py2:
        if not isinstance(resource, basestring):
            return resource
    else:
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


def create_random_string(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def gen_file_name(name=None, prefix='', output_dir='.', extension='log', append_time=True):
    """
    Generate a file name.

    :param name: a base for the file name
    :param output_dir: the directory
    :param extension: give the file an extension
    :param append_time: append time to the file name (before the extension)
    :return: the absolute path to the file
    """
    if name is None:
        name = create_random_string(15)

    if extension and extension[-1] == '.':
        extension = extension[-1:]

    if append_time:
        t = time.strftime("%Y-%m-%d.%H_%M_%S")
        file_name = "{}{}.{}.{}".format(prefix, name, t, extension)
    else:
        if extension:
            file_name = "{}{}.{}".format(prefix, name, extension)
        else:
            file_name = "{}{}".format(prefix, name)

    return os.path.abspath(os.path.join(output_dir, file_name))


def get_extension(filename):
    file_first, file_extension = os.path.splitext(filename)
    return file_extension


def prepend_before_extension(filename, text):
    file_first, file_extension = os.path.splitext(filename)

    if file_extension:
        return "{0}.{1}{2}".format(file_first, text, file_extension)

    return filename


def get_dir_and_file(filename):
    abspath = os.path.abspath(filename)
    if os.path.isdir(abspath):
        return abspath, None

    return os.path.split(abspath)


def check_file(filename, mode='r'):
    if mode == 'r':

        if filename and os.path.exists(filename):
            return os.path.abspath(filename)

        raise G2GValueError("The following file does not exist: {0}".format(filename))

    elif mode == 'w':
        file_dir = '.'

        if filename:
            file_name = os.path.abspath(filename)
            file_dir = os.path.dirname(file_name)

            if not os.access(file_dir, os.W_OK | os.X_OK):
                raise G2GValueError("Cannot generate file: {0}".format(filename))

            return file_name

    raise G2GValueError("Unspecified mode to open file, '{}'".format(mode))


def create_temp_dir(name, prefix='.g2gtools_', dir=None):
    new_name = "{}{}".format(prefix, name)
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
            raise G2GValueError("Will not delete directory: {}".format(dir))
        try:
            shutil.rmtree(dir)
        except Exception as e:
            raise G2GValueError("Will not delete directory: {}".format(dir))


def location_to_filestring(location):
    return "{}-{}-{}".format(location.seq_id, location.start, location.end)

#######################



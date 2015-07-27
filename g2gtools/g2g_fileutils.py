# -*- coding: utf-8 -*-

#
# Collection of g2g utility functions and classes
#

import bz2
import gzip
import os
import time
import urllib
from subprocess import Popen

from pysam import tabix_compress, faidx, tabix_index

from . import G2GValueError

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


def bgzip_index(original_file, new_file, file_format):
    """

    :param original_file:
    :param new_file:
    :param file_format:
    :return:
    """

    if file_format.lower() == 'fa':
        tabix_compress(original_file, new_file)
        faidx(new_file)
        delete_file(original_file)
    elif file_format.lower() == 'vcf':
        tabix_index(original_file, preset="vcf", force=True)
    else:
        raise G2GValueError("Unknown file format: {0}".format(file_format))


def open_resource(resource, mode='rb'):
    """
    Open different types of files and return the handle.

    :param resource: a file located locally or on the internet.  Gzip'd and zip'd files are handled.
    :param mode: standard file open modes
    :return: the resource (file) handle
    """
    if not resource:
        return None

    if not isinstance(resource, basestring):
        return resource

    if resource.endswith(('.gz', '.Z', '.z')):
        return gzip.open(resource, mode)
    elif resource.endswith(('.bz', '.bz2', '.bzip2')):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(('http://', 'https://', 'ftp://')):
        return urllib.urlopen(resource)
    else:
        return open(resource, mode)


def gen_file_name(name, output_dir='.', extension='log', append_time=True):
    """
    Generate a file name.

    :param name: a base for the file name
    :param output_dir: the directory
    :param extension: give the fils le an extension
    :param append_time: append time to the file name (before the extension)
    :return: the absolute path to the file
    """
    if extension[-1] == '.':
        extension = extension[-1:]

    if append_time:
        t = time.strftime("%Y-%m-%d.%H_%M_%S")
        file_name = "{0}.{1}.{2}".format(name, t, extension)
    else:
        file_name = "{0}.{1}".format(name, extension)

    return os.path.abspath(os.path.join(output_dir, file_name))



def get_extension(filename):
    file_first, file_extension = os.path.splitext(filename)
    return file_extension


def prepend_before_extension(filename, text):
    file_first, file_extension = os.path.splitext(filename)

    if file_extension:
        return "{0}.{1}{2}".format(file_first, text, file_extension)

    return filename


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


    raise G2GValueError("Unspecified mode")






if __name__ == '__main__':
    print ''

# -*- coding: utf-8 -*-

#
# Collection of functions related to Seqnature offset files
#

import os
import glob

from .g2g_utils import get_logger

LOG = get_logger()


def offset_parse_chromosomes(from_file, to_file):
    """
    Parse the chromosome information files and return a dictionary of chromosomes that contain information about the
    chromosome file offsets.

    It assumed that the chromosome file offsets are names "*_offsets_chr{chr}.txt"

    File format is one line per chromosome, with 2 fields: chromosome chromosome_length

    :param from_file: the chromosome file
    :param to_file: the chromosome file
    :return: a dictionary of chromosomes
    """
    chromosomes = {}

    try:
        if not os.path.exists(from_file):
            raise IOError(from_file)

        fd = open(from_file, 'r')
        for line in fd:
            elem = line.strip().split()
            chromosome = elem[0]
            chromosome_length = int(elem[1])
            chromosome_file_path = os.path.dirname(os.path.abspath(from_file))

            files = glob.glob1(chromosome_file_path, "*_offsets_chr{0}.txt".format(chromosome))
            chromosome_file_path = os.path.join(chromosome_file_path, files[0]) if len(files) == 1 else None
            chromosomes[chromosome] = {'from_chr':chromosome, 'from_length':chromosome_length, 'file_path':chromosome_file_path}
        fd.close()

        if not os.path.exists(to_file):
            raise IOError(to_file)

        fd = open(to_file, 'r')
        for line in fd:
            elem = line.strip().split()
            chromosome = elem[0]

            if chromosome in chromosomes:
                chromosome_length = int(elem[1])
                chromosomes[chromosome]['to_chr'] = chromosome
                chromosomes[chromosome]['to_length'] = chromosome_length
        fd.close()
    except IOError, e:
        message = "Error parsing chromosome files: {0}".format(e.message)
        LOG.info(message)
        return {}

    return chromosomes






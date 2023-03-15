#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_g2gtools
----------------------------------

Tests for `g2gtools` module.
"""

import unittest
# import os
# import shutil
# import pysam
# from . import bamsam, chain


class TestG2gtools(unittest.TestCase):

    def setUp(self):
        # self.path = os.path.realpath(__file__)
        # self.dir = os.path.dirname(self.path)
        # self.dir_parent = os.path.abspath(os.path.join(self.dir, os.pardir))
        #
        # self.chain_file = os.path.join(self.dir_parent, "example/example1.chain")
        # self.input_file = os.path.join(self.dir_parent, "example/example1.sam")
        # self.output_dir = os.path.join(self.dir_parent, "test_example1")
        #
        # os.mkdir(self.output_dir)
        # self.bed_file_out = os.path.join(self.output_dir, "example1.out.sam")
        pass

    def test_example1(self):
        """
        TODO: Not working yet

        :return:
        """
        # chain_file = chain.ChainFile(self.chain_file)
        # try:
        #     chain_file.parse()
        # except chain.ChainFileException, cfe:
        #     print cfe
        #
        # bamsam.convert_bam_file(chain=chain_file, infile=self.input_file, outfile=self.bed_file_out)
        #
        # try:
        #     sam_file_out = pysam.Samfile(self.bed_file_out, 'rb')
        # except:
        #     sam_file_out = pysam.Samfile(self.bed_file_out, 'r')
        #
        # try:
        #     while 1:
        #         old_alignment = sam_file_out.next()
        #         print old_alignment
        # except StopIteration:
        #     print("All reads processed")
        #
        # shutil.rmtree(self.output_dir)
        pass


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())

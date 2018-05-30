===============================
g2gtools
===============================

.. image:: https://badge.fury.io/py/g2gtools.png
    :target: http://badge.fury.io/py/g2gtools

.. image:: https://travis-ci.org/churchill-lab/g2gtools.png?branch=master
    :target: https://travis-ci.org/churchill-lab/g2gtools

.. image:: https://anaconda.org/kbchoi/g2gtools/badges/version.svg
    :target: https://anaconda.org/kbchoi/g2gtools

.. image:: https://readthedocs.org/projects/g2gtools/badge/?version=latest
    :target: http://g2gtools.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://zenodo.org/badge/39776319.svg
    :target: https://zenodo.org/badge/latestdoi/39776319

**g2gtools** creates custom genomes by incorporating SNPs and indels into reference genome, extracts regions of interest, e.g., exons or transcripts, , and converts coordinates of files (bam, gtf, bed) between two genomes. Unlike other liftover tools, g2gtools does not throw away alignments that land on indel regions. Release Version 0.2 can now create personal **diploid** genomes. The new version still lifts over diploid alignments on personal genome coordinates back to that of reference so we can compare alignments from among samples in a population.

* Free software: MIT License
* Documentation for Release Ver. 0.2.XX: http://churchill-lab.github.io/g2gtools/.
* Documentation for Release Ver. 0.1.XX: https://g2gtools.readthedocs.org.


Reference
~~~~~~~~~

Manuscript in preparation (expected in 2018)
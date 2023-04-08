===============================
g2gtools
===============================

**g2gtools** creates custom genomes by incorporating SNPs and indels into reference genome, extracts regions of interest, e.g., exons or transcripts, , and converts coordinates of files (bam, gtf, bed) between two genomes. Unlike other liftover tools, g2gtools does not throw away alignments that land on indel regions. Release Version 0.2 can now create personal **diploid** genomes. The new version still lifts over diploid alignments on personal genome coordinates back to that of reference so we can compare alignments from among samples in a population.

* Free software: MIT License
* Documentation for Release Ver. 1.0.XX: BELOW
* Documentation for Release Ver. 0.2.XX: http://churchill-lab.github.io/g2gtools/.
* Documentation for Release Ver. 0.1.XX: https://g2gtools.readthedocs.org.

Development Lead
----------------

* Algorithm design and software engineering: **Kwangbom "KB" Choi, Ph. D.**, The Jackson Laboratory <kb.choi@jax.org>
* Software engineering: **Matthew J. Vincent**, The Jackson Laboratory <matt.vincent@jax.org>

Contributors
------------

* Narayanan Raghupathy, The Jackson Laboratory <Narayanan.Raghupathy@jax.org>
* Anuj Srivastava, The Jackson Laboratory <Anuj.Srivastava@jax.org>
* Mike Lloyd, The Jackson Laboratory <Mike.Lloyd@jax.org>

Docker
------------
A docker file is provided for easy build.  We will also be pushing to DockerHub.

The main portion of this release was to upgrade to Python 3.x.  Please reach out with any bugs you may encounter.


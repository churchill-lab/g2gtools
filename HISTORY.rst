.. :changelog:

History
-------

3.0.0 (02/16/2025)
~~~~~~~~~~~~~~~~~~

* Revamped all the commands and logging to be more consistent and user-friendly.
* Updated third party libraries.
* Tested on Python 3.10+.

2.4.0 (02/15/2025)
~~~~~~~~~~~~~~~~~~

* Fixed a bug in parsing GTF file where the attribute field's values contained nested semi-colons.

2.3.0 (02/03/2025)
~~~~~~~~~~~~~~~~~~

* Fixed a bug in parsing regions that have a period in their identifier.

2.2.0 (01/29/2025)
~~~~~~~~~~~~~~~~~~

* Minor fixes and some bug fix for GitHb issue #32


2.0.0 (04/06/2023)
~~~~~~~~~~~~~~~~~~

* Upgraded to Python 3 (dropped Python 2)

0.2.9 (10/01/2019)
~~~~~~~~~~~~~~~~~~

* Fixed error of not setting correct header information in B/SAM.

0.2.8 (09/23/2019)
~~~~~~~~~~~~~~~~~~

* Fixed error of not parsing VCI file on bam/sam conversion.

0.2.5 (06/06/2018)
~~~~~~~~~~~~~~~~~~

* Fixed error of always generating index file

0.2.4 (06/05/2018)
~~~~~~~~~~~~~~~~~~

* Automtically generates file index if not found

0.2.3 (06/05/2018)
~~~~~~~~~~~~~~~~~~

* Fixed end of contig tree mapping

0.2.2 (06/04/2018)
~~~~~~~~~~~~~~~~~~

* Fixed extract transcripts

0.2.1 (05/30/2018)
~~~~~~~~~~~~~~~~~~

* Hotfix to solve gtf parsing problem in python 3

0.2.0 (05/09/2018)
~~~~~~~~~~~~~~~~~~

* Support for diploid VCF
* Reimplemented to support G2G VCI (Variant Call Information)
* Support python3

0.1.31 (02/17/2017)
~~~~~~~~~~~~~~~~~~~

* Final release for custom haploid reconstruction

0.1.29 (05/17/2016)
~~~~~~~~~~~~~~~~~~~

* Fixed parsing issue with GATK-generated VCF files

0.1.27 (05/04/2016)
~~~~~~~~~~~~~~~~~~~

* Uploaded the package to Anaconda Cloud
* Fixed Travis CI fail

0.1.23 (02/15/2016)
~~~~~~~~~~~~~~~~~~~

* Setup Travis CI for automated building

0.1.22 (02/15/2016)
~~~~~~~~~~~~~~~~~~~

* Updated documentation

0.1.20 (02/03/2016)
~~~~~~~~~~~~~~~~~~~

* Released to public with documentation at g2gtools.readthedocs.org

0.1.0 (02/09/2015)
~~~~~~~~~~~~~~~~~~

* Started github repository

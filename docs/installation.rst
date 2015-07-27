============
Installation
============

We recommend using conda environment::

    $ git clone https://github.com/churchill-lab/g2gtools.git
    $ conda create -n g2gtools pip biopython>=1.63
    $ source activate g2gtools
    (g2gtools) $ conda install -c https://conda.binstar.org/bcbio pysam
    (g2gtools) $ conda install -c https://conda.binstar.org/bcbio bx-python
    (g2gtools) $ python setup.py install

If you have virtualenvwrapper installed::

    $ mkvirtualenv g2gtools
    $ pip install g2gtools

Or, at the command line::

    $ easy_install g2gtools


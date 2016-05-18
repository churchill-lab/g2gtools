============
Installation
============

Although g2gtools is available at PyPI (https://pypi.python.org/pypi/g2gtools/) for 'pip install' or 'easy_install', we highly recommend using Anaconda distribution of python (https://www.continuum.io/downloads) to install all the dependencies without issues. The most recent version of g2gtools is also available on Anaconda Cloud, so add the following channels if you have not already::

    $ conda config --add channels r
    $ conda config --add channels bioconda

To avoid conflicts among dependencies, we highly recommend using conda virtual environment::

    $ conda create -n g2gtools jupyter
    $ source activate g2gtools

Once g2gtools virtual environment is created and activated, your shell prompt will show '(g2gtools)' at the beginning to specify what virtual environment you are currently in. Now please type the following and install g2gtools::

    (g2gtools) $ conda install -c kbchoi g2gtools

That's all! We note that you can go out from g2gtools virtual environment anytime once you are done using g2gtools::

    (g2gtools) $ source deactivate


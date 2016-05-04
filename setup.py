#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.extension import Extension

from setuptools import *
from glob import glob


try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

from distutils.core import Command

# Use build_ext from Cython
command_classes = {}

# Use build_ext from Cython if found
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
except:
    pass

# Run 2to3 builder if we're on Python 3.x, from
#   http://wiki.python.org/moin/PortingPythonToPy3k
try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    # 2.x
    from distutils.command.build_py import build_py
command_classes['build_py'] = build_py

requirements = []

#with open('requirements.txt') as requirements_file:
#    for line in requirements_file:
#        requirements.append(line.strip())

test_requirements = [
    'pytest'
]


import os
on_rtd = os.environ.get('READTHEDOCS', None)

if not on_rtd:
    requirements.append('bx-python>=0.7.2')
    requirements.append('pysam>=0.8.3')
    requirements.append('biopython>=1.63')
    requirements.append('pysqlite>=2.6.3')

import g2gtools


def main():
    with open('README.rst') as readme_file:
        readme = readme_file.read()

    with open('HISTORY.rst') as history_file:
        history = history_file.read().replace('.. :changelog:', '')

    setup(
        name='g2gtools',
        version=g2gtools.__version__,
        description="A set of tools that facilitates genome to genome conversion for studying multiparent populations",
        long_description=readme + '\n\n' + history,
        author='Kwangbom "KB" Choi & Matthew Vincent, The Jackson Laboratory',
        author_email='kb.choi@jax.org',
        url='https://github.com/churchill-lab/g2gtools',
        packages=find_packages('.'),
        package_dir={'g2gtools': 'g2gtools'},
        scripts=glob("bin/*"),
        ext_modules=get_extension_modules(),
        cmdclass=command_classes,
        include_package_data=True,
        install_requires=requirements,
        license="GPLv3",
        keywords='g2gtools',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Programming Language :: Python :: 2.7'
        ],
        test_suite='tests',
        tests_require=test_requirements,
        zip_safe=False
    )


def get_extension_modules():
    extensions = []

    # fasta reading speedups
    #extensions.append(Extension("g2gtools._fasta", ["lib/g2gtools/_fasta.pyx", "lib/g2gtools/Fasta.cpp"], language='c++'))

    return extensions

if __name__ == '__main__':
    main()

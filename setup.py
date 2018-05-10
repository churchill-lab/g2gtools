#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from glob import glob
with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

import os
on_rtd = os.environ.get('READTHEDOCS', None)

requirements = []
test_requirements = []

if not on_rtd:
    requirements.append('Cython')
    requirements.append('numpy>=1.14')
    requirements.append('bx-python>=0.8')
    requirements.append('pysam>=0.14')
    requirements.append('natsort>=5.0.1')
    requirements.append('future>=0.15')

setup(
    name='g2gtools',
    version='0.2.0',
    description="A set of tools that facilitate the reconstruction of custom diploid genomes and coordinate conversion",
    long_description=readme + '\n\n' + history,
    author='Kwangbom "KB" Choi & Matthew Vincent, The Jackson Laboratory',
    author_email='kb.choi@jax.org',
    url='https://github.com/churchill-lab/g2gtools',
    packages=[
        'g2gtools',
    ],
    package_dir={'g2gtools':
                 'g2gtools'},
    scripts=glob("bin/*"),
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3",
    zip_safe=False,
    keywords='g2gtools',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7'
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.5'
        'Programming Language :: Python :: 3.6'
    ],
    test_suite='tests',
    tests_require=test_requirements
)

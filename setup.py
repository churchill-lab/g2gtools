#!/usr/bin/env python
import os
from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

requirements = []
test_requirements = []

on_rtd = os.environ.get("READTHEDOCS", None)

if not on_rtd:
    with open("requirements.txt") as requirements_file:
        requirements_lines = requirements_file.readlines()
        for line in requirements_lines:
            requirements.append(line)

setup(
    name="g2gtools",
    version="2.0.0",
    description="A suite of tools for the reconstruction of personal diploid genomes and better coordinate conversion",
    long_description=readme,
    author="Matthew J. Vincent and Kwangbom 'KB' Choi, The Jackson Laboratory",
    author_email="matt.vincent@jax.org",
    url="https://github.com/churchill-lab/g2gtools",
    keywords=["g2gtools", "personal genomes", "liftover"],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "g2gtools = g2gtools.commands:cli",
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.10",
    ],
)

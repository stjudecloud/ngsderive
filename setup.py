#!/usr/bin/env python

from setuptools import setup, find_packages

__VERSION__ = "1.0.1"

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open("README.md", "r") as fh:
    readme = fh.read()

setup(
  name='ngsderive',
  version=__VERSION__,
  description='Backwards derive attributes from NGS data',
  long_description=readme,
  long_description_content_type="text/markdown",
  author='Clay McLeod',
  author_email='clay.mcleod@STJUDE.org',
  scripts=['scripts/ngsderive'],
  packages=find_packages(),
  install_requires=requirements,
  python_requires='>=3.0, <3.8',
  classifiers = [
      "Programming Language :: Python :: 3",
      "Operating System :: OS Independent",
      "Development Status :: 5 - Production/Stable",
      "Topic :: Scientific/Engineering :: Bio-Informatics"
  ]
)

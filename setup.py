#!/usr/bin/env python

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
  name='ngsderive',
  version='1.0.0',
  description='Backwards derive attributes from NGS data',
  author='Clay McLeod',
  author_email='clay.mcleod@STJUDE.org',
  scripts=['scripts/ngsderive'],
  packages=find_packages(),
  install_requires=requirements,
  python_requires='>=3.0, <=3.7'
)

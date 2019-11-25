#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
  name='ngsderive',
  version='1.0.0',
  description='Backwards derive information from NGS data',
  author='Clay McLeod',
  author_email='clay.mcleod@STJUDE.org',
  scripts=['scripts/ngsderive'],
  packages=find_packages(),
)

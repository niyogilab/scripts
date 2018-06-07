#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'qrlabels',
  version = '0.1',
  packages = ['qrlabels'],
  entry_points = {'console_scripts': ['qrlabels=qrlabels:main']}
)

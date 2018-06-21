#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'tecan_extract_table',
  version = '1.0.0',
  packages = ['tecan_extract_table'],
  entry_points = {'console_scripts': ['tecan-extract-table=tecan_extract_table:main']},
  package_data = {'tecan_extract_table': ['usage.txt']}
)

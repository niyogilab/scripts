#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'tidytecan',
  version = '0.1.1',
  packages = ['tidytecan'],
  entry_points = {'console_scripts': ['tidytecan=tidytecan:main']},
  package_data = {'tidytecan': ['usage.txt']}
)

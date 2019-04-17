#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'cyano_cpf1_primers',
  version = '1.0',
  packages = ['cyano_cpf1_primers'],
  entry_points = {'console_scripts': ['cyano_cpf1_primers=cyano_cpf1_primers:main']}
)

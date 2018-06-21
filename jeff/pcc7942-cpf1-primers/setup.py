#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'pcc7942_cpf1_primers',
  version = '1.0',
  packages = ['pcc7942_cpf1_primers'],
  entry_points = {'console_scripts': ['pcc7942-cpf1-primers=pcc7942_cpf1_primers:main']}
)

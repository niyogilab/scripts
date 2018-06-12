#!/usr/bin/env python

from setuptools import setup

setup(
  name = 'cpf1primers',
  version = '0.1',
  packages = ['cpf1primers'],
  entry_points = {'console_scripts': ['cpf1primers=cpf1primers:main']}
)

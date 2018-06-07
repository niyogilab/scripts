#!/usr/bin/env python

'''
Generates a table of initial primers to try for Cpf1 knockouts of PCC 7942 genes.

Usage:
  cpf1primers (-h | --help)
  cpf1primers [-v] LOCUS...

Options:
  -h, --help  Show this text
  -v          Print debugging information to stderr
  LOCUS       Locus ID to generate primers for. You can put more than one.
              For example, 0001 1021 means Synpcc7942_0001 and Synpcc7942_1021.
'''

from __future__ import print_function
from docopt     import docopt
from sys        import stderr

### interface ###

def log(msg, args):
  if args['verbose']:
      print(args, file=stderr)

def parse(args):
  return {
    'verbose'  : args['-v'],
    'locusids' : ['Synpcc7942_%s' % a for a in args['LOCUS']]
  }

def main():
  args = parse(docopt(__doc__, version='cpf1primers 0.1'))
  log(args, args)

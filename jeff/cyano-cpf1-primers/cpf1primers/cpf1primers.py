#!/usr/bin/env python

# TODO fix errors at beginning/end of chr, and allow pANL genes?

'''
Generates a table of initial primers to try for Cpf1 knockouts of PCC 7942 genes.

Usage:
  cpf1primers (-h | --help)
  cpf1primers [-v] [-g GENOME] LOCUS...

Options:
  -h, --help  Show this text
  -v          Print debugging information to stderr
  -g GENOME   PCC 7942 genome to use. [default: Synpcc7942_chr.gbk]
  LOCUS       Locus ID to generate primers for. You can put more than one.
              For example, 0001 1021 means Synpcc7942_0001 and Synpcc7942_1021.
'''

from __future__ import print_function
from docopt     import docopt
from sys        import stderr

def get_coordinates(args, locus):
    coords = (1000, 2000) # TODO write this
    log(args, '%s coordinates: %s' % (locus, coords))
    return coords

def get_sequence(args, locus):
    sequence = 'acagacgtactgatcgtagctacgtacgtacggtcgcc' # TODO write this
    log(args, '%s sequence: %s' % (locus, sequence))
    return sequence

def find_targets(args, locus, sequence):
    # TODO also need to find reverse targets?
    targets = ['tggatcgatcgatcgatgcgatcgtacgtacg']
    log(args, 'found %s potential targets in %s: %s' % (len(targets), locus, targets))
    return targets

def guide_rna(args, locus):
    sequence = get_sequence(args, locus)
    pam_start = find_targets(args, locus, sequence)
    log(args, '%s guide rna   forward primer' % locus)
    log(args, '%s guide rna   reverse primer' % locus)

def hr_template(args, locus, homology_bp=750):
    gene_start , gene_end  = get_coordinates(args, locus)
    left_start  = gene_start - homology_bp
    left_end    = gene_start
    right_start = gene_end
    right_end   = gene_end + homology_bp
    log(args, '%s left  flank (%s-%sbp) forward primer' % (locus, left_start , left_end))
    log(args, '%s left  flank (%s-%sbp) reverse primer' % (locus, left_start , left_end))
    log(args, '%s right flank (%s-%sbp) forward primer' % (locus, right_start, right_end))
    log(args, '%s right flank (%s-%sbp) reverse primer' % (locus, right_start, right_end))

def log(args, msg):
  if args['verbose']:
      print(msg, file=stderr)

def parse(args):
  return {
    'verbose'  : args['-v'],
    'genome'   : args['-g'],
    'locusids' : ['Synpcc7942_%s' % a for a in args['LOCUS']]
  }

def main():
  args = parse(docopt(__doc__, version='cpf1primers 0.1'))
  log(args, args)
  for locus in args['locusids']:
      guide_rna(args, locus)
      hr_template(args, locus)

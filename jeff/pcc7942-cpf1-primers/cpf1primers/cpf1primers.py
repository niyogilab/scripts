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
  -g GENOME   PCC 7942 genome to use. [default: SynPCC7942_chr.gbk]
  LOCUS       Locus ID to generate primers for. You can put more than one.
              For example, 0001 1021 means Synpcc7942_0001 and Synpcc7942_1021.
'''

from __future__ import print_function
from Bio        import SeqIO
from Bio.Seq    import Seq
from Bio.Alphabet import generic_dna
from docopt     import docopt
from sys        import stderr
import re

# silence warnings
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

def seqid(seq_feature):
  for q in ['locus_tag', 'protein_id']:
    if q in seq_feature.qualifiers:
      return seq_feature.qualifiers[q][0] # TODO is there always 1?
  raise Exception('no seqid found for feature: %s' % seq_feature)

def get_coordinates(args, locus):
    coords = (1000, 2000) # TODO write this
    log(args, '%s coordinates: %s' % (locus, coords))
    return coords

def get_sequences(args):
  loci = args['locusids']
  seqs = []
  for seq_record in SeqIO.parse(args['genome'], 'genbank'):
    for seq_feature in seq_record.features:
      if seq_feature.type != 'CDS':
        continue
    # adapted from:
    # http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta
    # for seq_record in SeqIO.parse(ingbk, 'genbank'):
      # for seq_feature in seq_record.features:
        # if seq_feature.type != 'CDS': # TODO anything else we might want?
          # continue
        # assert len(seq_feature.qualifiers['locus_tag']) == 1
        # assert len(seq_feature.qualifiers['product']) == 1
        # assert len(seq_feature.qualifiers['translation']) == 1
      sid = seqid(seq_feature)
      #log(args, 'checking seqid %s' % sid)
      if sid in loci:
        seqstr = str(seq_feature.extract(seq_record.seq))
        seqs.append({'locusid': sid, 'sequence': seqstr})

    sequence = 'acagacgtactgatcgtagctacgtacgtacggtcgcc' # TODO write this
    log(args, 'sequences found: %s' % seqs)
    # log(args, '%s sequence: %s' % (loci, sequence))
    return seqs

def find_targets(args, locus, sequence):
  # TODO also need to find reverse targets?
  pam_and_target = '[TC]T.{21}'
  targets = re.findall(pam_and_target, sequence)
  log(args, 'found %s potential targets in %s: %s' % (len(targets), locus, str(targets)))
  return targets

def guide_rna(args, seq):
    locus = seq['locusid']
    # sequence = get_sequence(args, locus)
    targets = find_targets(args, locus, seq['sequence'])
    target  = targets[0] # TODO any reason to be more selective?
    primer_fwd = 'AGAT' + target[2:]
    primer_rev = 'AGAC' + str(Seq(target[2:], generic_dna).reverse_complement())
    log(args, '%s guide rna forward primer: %s' % (locus, primer_fwd))
    log(args, '%s guide rna reverse primer: %s' % (locus, primer_rev))
    return (primer_fwd, primer_rev)

def hr_template(args, seq, homology_bp=750):
  locus = seq['locusid']
  gene_start , gene_end = get_coordinates(args, locus)
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
  seqs = get_sequences(args)
  headers = ['locus', 'guide_fwd', 'guide_rev']
  print('\t'.join(headers))
  for seq in seqs:
      guide_fwd, guide_rev = guide_rna(args, seq)
      row = (seq['locusid'], guide_fwd, guide_rev)
      print('\t'.join(row))
      # hr_template(args, seq)

#!/usr/bin/env python

# TODO fix errors at beginning/end of chr, and allow pANL genes?

'''
Generates a table of initial primers to try for Cpf1 knockouts of PCC 7942 genes.

Usage:
  pcc7942-cpf1-primers (-h | --help)
  pcc7942-cpf1-primers [-v] [-g GENOME] LOCUS...

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

def get_sequences(args, genome):
  loci = args['locusids']
  seqs = []
  for seq_feature in genome.features:
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
      seqstr = str(seq_feature.extract(genome.seq))
      start  = int(seq_feature.location.start)
      end    = int(seq_feature.location.end)
      seqs.append({'locusid': sid, 'sequence': seqstr, 'start': start, 'end': end})

    # sequence = 'acagacgtactgatcgtagctacgtacgtacggtcgcc' # TODO write this
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

def hr_primers(args, genome, seq, homology_bp=750):
  locus = seq['locusid']

  # get flanking sequences
  left_start  = seq['start'] - homology_bp
  left_end    = seq['start']
  right_start = seq['end']
  right_end   = seq['end'] + homology_bp

  # design primers for them
  # note this is just a basic first pass! TODO use primer3 here
  left_fwd  = str(genome.seq[left_start  : left_start  + 20])
  right_fwd = str(genome.seq[right_start : right_start + 20].reverse_complement())
  left_rev  = str(genome.seq[left_end  - 20 : left_end ])
  right_rev = str(genome.seq[right_end - 20 : right_end].reverse_complement())

  log(args, '%s left  flank (%s-%sbp) forward primer: %s' % (locus, left_start , left_end, left_fwd))
  log(args, '%s left  flank (%s-%sbp) reverse primer: %s' % (locus, left_start , left_end, left_rev))
  log(args, '%s right flank (%s-%sbp) forward primer: %s' % (locus, right_start, right_end, right_fwd))
  log(args, '%s right flank (%s-%sbp) reverse primer: %s' % (locus, right_start, right_end, right_rev))

  return [left_fwd, left_rev, right_fwd, right_rev]

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
  genome = list(SeqIO.parse(args['genome'], 'genbank'))[0] # TODO look up the right way
  log(args, 'genom: %s' % str(genome.seq)[:1000])
  seqs = get_sequences(args, genome)
  headers = ['locus', 'guide_fwd', 'guide_rev', 'left_fwd', 'left_rev', 'right_fwd', 'right_rev']
  print('\t'.join(headers))
  for seq in seqs:
      row = [seq['locusid']]
      guide_fwd, guide_rev = guide_rna(args, seq)
      row += [guide_fwd, guide_rev]
      # left_fwd, left_rev, right_fwd, right_rev = hr_primers(args, seq)
      row += hr_primers(args, genome, seq)
      print('\t'.join(row))

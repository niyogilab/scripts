#!/usr/bin/env python3

'''
Generates a table of initial primers to try for Cpf1 knockouts of PCC 7942 genes.

Designed for markerless KOs with pSL2680 from the Pakrasi lab (AddGene plasmid #85581).
The crRNA pair should be annealed and ligated into the plasmid directly;
left and right pairs are for generating a standard HR template by PCR with the target gene missing.
The final HR template has KpnI and SalI sites for integrating into pSL2680.

Usage:
  pcc7942_cpf1_primers (-h | --help)
  pcc7942_cpf1_primers [-v] [-g GENOME] LOCUS...

Options:
  -h, --help  Show this text
  -v          Print debugging information to stderr
  -g GENOME   PCC 7942 genome to use. [default: SynPCC7942_chr.gbk]
  LOCUS       Locus ID to generate primers for. You can put more than one.
              For example, 0001 1021 means Synpcc7942_0001 and Synpcc7942_1021.
'''

from __future__  import print_function

import primer3
import re

from Bio          import SeqIO
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from copy         import deepcopy
from docopt       import docopt
from sys          import stderr

# silence warnings
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

def seqid(seq_feature):
  for q in ['locus_tag', 'protein_id']:
    if q in seq_feature.qualifiers:
      return seq_feature.qualifiers[q][0] # TODO is there always 1?
  raise Exception('no seqid found for feature: %s' % seq_feature)

# def get_coordinates(args, locus):
#     coords = (1000, 2000) # TODO write this
#     log(args, '%s coordinates: %s' % (locus, coords))
#     return coords

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
  # TODO should the targets actually include the end of the left flank so they start immediately?
  pam_and_target = '[TC]T.{21}'
  targets = re.findall(pam_and_target, sequence)
  log(args, 'found %s potential CPF1 targets in %s: %s' % (len(targets), locus, str(targets)))
  # TODO error if 0
  return targets

def guide_rna(args, seq):
    locus = seq['locusid']
    # sequence = get_sequence(args, locus)
    targets = find_targets(args, locus, seq['sequence'])
    # TODO any reason to be more selective?
    # TODO option to pick more than one per locus?
    # TODO can i use primer3 here to to guess which ones will be more effective?
    target  = targets[0]
    log(args, 'using the first one, %s' % target)
    primer_fwd = 'AGAT' + target[2:]
    primer_rev = 'AGAC' + str(Seq(target[2:], generic_dna).reverse_complement())
    log(args, '%s guide rna forward primer: %s' % (locus, primer_fwd))
    log(args, '%s guide rna reverse primer: %s' % (locus, primer_rev))
    return (primer_fwd, primer_rev)

# doing this once at the beginning is theoretically faster, but mostly simple
# see https://libnano.github.io/primer3-py/quickstart.html#primer-design
# TODO fiddle with these or find sane defaults on the web
PRIMER_ARGS = {

  # pick left + right primers (TODO different for guide rna?)
  'PRIMER_TASK': 'generic',
  'PRIMER_PICK_LEFT_PRIMER'   : True,
  'PRIMER_PICK_INTERNAL_OLIGO': False,
  'PRIMER_PICK_RIGHT_PRIMER'  : True,

  # 'PRIMER_OPT_SIZE': 20,
  # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
  # 'PRIMER_MIN_SIZE': 18,
  # 'PRIMER_MAX_SIZE': 25,
  # 'PRIMER_OPT_TM': 60.0,
  # 'PRIMER_MIN_TM': 57.0,
  # 'PRIMER_MAX_TM': 63.0,
  # 'PRIMER_MIN_GC': 20.0,
  # 'PRIMER_MAX_GC': 80.0,
  # 'PRIMER_MAX_POLY_X': 100,
  # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
  # 'PRIMER_SALT_MONOVALENT': 50.0,
  # 'PRIMER_DNA_CONC': 50.0,
  # 'PRIMER_MAX_NS_ACCEPTED': 0,
  # 'PRIMER_MAX_SELF_ANY': 12,
  # 'PRIMER_MAX_SELF_END': 8,
  # 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
  # 'PRIMER_PAIR_MAX_COMPL_END': 8,
  'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],
                               [150,175],[175,200],[200,225]],
}

# TODO should i have a hybridization probe?
# TODO extract and return just sequences first
# TODO extract Tm too, and add to the table
def run_primer3(args, seqid, seq, extra_seq_args={}):
  extra_primer_args = {
    # 'PRIMER_PRODUCT_SIZE_RANGE': [[100, 10000]]
  }
  seq_args = {
    'SEQUENCE_ID': seqid,
    'SEQUENCE_TEMPLATE': seq,
  }
  seq_args.update(extra_seq_args)
  primer_args = deepcopy(PRIMER_ARGS)
  primer_args.update(extra_primer_args)
  log(args, 'primer_args: %s' % primer_args)
  log(args, 'seq_args: %s' % seq_args)
  primers = primer3.designPrimers(seq_args, primer_args)
  return primers

def primer3_left(args, seqid, seq):
  return run_primer3(args, seqid, seq, {'SEQUENCE_FORCE_LEFT_END': len(seq)}) # TODO one too many?

def primer3_right(args, seqid, seq):
  return run_primer3(args, seqid, seq, {'SEQUENCE_FORCE_RIGHT_START': 1})

# TODO primer3_grna

def extract_first_pair(args, p3res):
  'extract just relevant info about the first primer pair'
  return (p3res['PRIMER_LEFT_1_SEQUENCE' ],
          '%0.1f' % p3res['PRIMER_LEFT_1_TM' ],
          p3res['PRIMER_RIGHT_1_SEQUENCE'],
          '%0.1f' % p3res['PRIMER_RIGHT_1_TM'])

def hr_primers(args, genome, seq, homology_bp=1000):
  locus = seq['locusid']

  # get flanking sequences
  left_start  = seq['start'] - homology_bp
  left_end    = seq['start']
  right_start = seq['end']
  right_end   = seq['end'] + homology_bp

  # design primers for them
  left_res  = run_primer3(args, 'left' , str(genome.seq[left_start  :  left_end]))
  right_res = run_primer3(args, 'right', str(genome.seq[right_start : right_end]))
  left_fwd_seq , left_fwd_tm , left_rev_seq , left_rev_tm  = extract_first_pair(args, left_res )
  right_fwd_seq, right_fwd_tm, right_rev_seq, right_rev_tm = extract_first_pair(args, right_res)

  log(args, '%s left  flank (%s-%sbp) forward primer: %s (Tm=%s)' % (locus, left_start , left_end , left_fwd_seq , left_fwd_tm ))
  log(args, '%s left  flank (%s-%sbp) reverse primer: %s (Tm=%s)' % (locus, left_start , left_end , left_rev_seq , left_rev_tm ))
  log(args, '%s right flank (%s-%sbp) forward primer: %s (Tm=%s)' % (locus, right_start, right_end, right_fwd_seq, right_fwd_tm))
  log(args, '%s right flank (%s-%sbp) reverse primer: %s (Tm=%s)' % (locus, right_start, right_end, right_rev_seq, right_rev_tm))

  return [left_fwd_seq,  left_fwd_tm,  left_rev_seq,  left_rev_tm,
         right_fwd_seq, right_fwd_tm, right_rev_seq, right_rev_tm]

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
  # log(args, 'genome: %s' % str(genome.seq)[:1000])
  seqs = get_sequences(args, genome)
  headers = ['locus_ID', 'crRNA_fwd', 'crRNA_rev',
              'left_fwd_seq',  'left_fwd_tm', ' left_rev_seq',  'left_rev_tm',
             'right_fwd_seq', 'right_fwd_tm', 'right_rev_seq', 'right_rev_tm']
  print('\t'.join(headers))
  for seq in seqs:
      row = [seq['locusid']]
      guide_fwd, guide_rev = guide_rna(args, seq)
      row += [guide_fwd, guide_rev]
      # left_fwd, left_rev, right_fwd, right_rev = hr_primers(args, genome, seq)
      row += hr_primers(args, genome, seq)
      print('\t'.join(row))

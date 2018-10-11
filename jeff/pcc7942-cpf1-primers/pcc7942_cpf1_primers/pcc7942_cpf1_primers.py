#!/usr/bin/env python3

# TODO refactor to pass around dicts rather than tuples
# TODO refactor to print the table in IDT order format
# TODO log levels

'''
Generates a table of initial primers to try for Cpf1 knockouts of PCC 7942 genes.

Designed for markerless KOs with pSL2680 from the Pakrasi lab (AddGene plasmid #85581).
The crRNA pair should be annealed and ligated into the plasmid directly;
left and right pairs are for generating a standard HR template by PCR with the target gene missing.
The final HR template has KpnI and SalI sites for integrating into pSL2680.

Usage:
  pcc7942_cpf1_primers (-h | --help)
  pcc7942_cpf1_primers [-v...] [-g GENOME] LOCUS...

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

def log(args, msg, level=1):
  if args['verbose'] >= level:
    print(msg, file=stderr)

def shorten_seqs(args, seq_dict):
  'this is just to make printing them less horrible'
  shorter = {}
  for key in seq_dict:
    val = str(seq_dict[key])
    if len(val) > 20:
      val = val[:20] + '...'
    shorter[key] = val
  return shorter


################################
# find sequences in the genome #
################################

def find_seqid(seq_feature):
  'takes a sequence feature and returns its locus id. exception if none'
  for q in ['locus_tag', 'protein_id']:
    if q in seq_feature.qualifiers:
      return seq_feature.qualifiers[q][0] # TODO is there always 1?
  raise Exception('no seqid found for feature: %s' % seq_feature)

def find_sequences(args, genome, loci, homology_bp=1000):
  '''
  takes a genome and a list of loci.
  returns a list of dicts with a left, right and cds sequence for each locus.
  '''
  # http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta
  seqs = []
  for seq_feature in genome.features:
    # skip if not one of the requested loci
    if seq_feature.type != 'CDS':
      continue
    sid = find_seqid(seq_feature)
    if not sid in loci:
      continue
    # define the gene itself
    start = int(seq_feature.location.start)
    end   = int(seq_feature.location.end)
    # define the flanking homology arms
    left_start  = start - homology_bp
    left_end    = start
    right_start = end
    right_end   = end + homology_bp
    # TODO is orientation important here, or can we CRISPR it in reverse just as well?
    seqs.append({
      'id'    : sid,
      'cds'   : str(seq_feature.extract(genome.seq)), # TODO same as indexing start : end right?
      'left'  : str(genome.seq[left_start  :  left_end]),
      'right' : str(genome.seq[right_start : right_end]),
    })
  log(args, 'seqs: %s' % [shorten_seqs(args, s) for s in seqs], 2)
  return seqs

##############################
# generate pre-crRNA primers #
##############################

# TODO no need to find Tm if being at the beginning is more important?
# TODO but should look for potential off-target sites

def find_cpf1_targets(args, locus, sequence):
  # TODO also need to find reverse targets?
  # TODO should the targets actually include the end of the left flank so they start immediately?
  pam_and_target = '[TC]T.{21}'
  targets = [t[2:] for t in re.findall(pam_and_target, sequence)]
  log(args, 'found %s potential CPF1 targets in %s: %s' % (len(targets), locus, str(targets)), 2)
  # TODO error if 0
  return targets

# this works, but i think having the cut site next to the flank might be more important
# def pick_lowest_tm(args, pairs):
#   'given lists of (fwd, rev) primers pairs, return the one with the lowest Tm, plus its Tm'
#   low_pick = (None, None, 999)
#   for (fwd, rev) in pairs:
#     tm = primer3.calcHeterodimerTm(fwd, rev)
#     if tm < low_pick[2]:
#       low_pick = (fwd, rev, tm)
#   return low_pick

def crRNA_primers_for_cpf1_target(target):
  'generates fwd and rev primers with AarI sites added'
  return {
    'fwd': 'AGAT' + target,
    'rev': 'AGAC' + str(Seq(target, generic_dna).reverse_complement())
  }

def crRNA(args, seq):
  locus = seq['id']
  # sequence = get_sequence(args, locus)
  targets = find_cpf1_targets(args, locus, seq['cds'])
  primers = crRNA_primers_for_cpf1_target(targets[0])

  # TODO any reason to be more selective?
  # TODO option to pick more than one per locus?
  # TODO can i use primer3 here to to guess which ones will be more effective?
  # target  = targets[0]
  log(args, 'using the first one: %s' % primers, 2)
  # fwd_primers = ['AGAT' + t for t in targets]
  # rev_primers = ['AGAC' + str(Seq(t, generic_dna).reverse_complement()) for t in targets]
  # pairs = zip(fwd_primers, rev_primers)

  # this finds the lowest-Tm pair, but I'm not sure if being nearest the flank is more important
  #(fwd_primer, rev_primer, tm) = pick_lowest_tm(args, pairs))
  # log(args, '%s guide rna annealing temp: %s' % (locus, tm))
  # (fwd_primer, rev_primer) = pairs[0]
  # fwd_primer = primers['fwd']
  # rev_primer = primers['rev']

  log(args, '%s crRNA forward primer: %s' % (locus, primers['fwd']), 2)
  log(args, '%s crRNA reverse primer: %s' % (locus, primers['rev']), 2)
  return primers


#######################################
# generate left + right flank primers #
#######################################

# TODO fiddle with these or find sane defaults on the web
# TODO should i have a hybridization probe?
# TODO extract and return just sequences first
# TODO extract Tm too, and add to the table
def suggest_primers(args, seq, extra_seq_args={}):
  seq_args = {
    'SEQUENCE_ID': 'seqid', # TODO any reason to have use this?
    'SEQUENCE_TEMPLATE': seq,
  }
  seq_args.update(extra_seq_args)
  primer_args = {
    # see https://libnano.github.io/primer3-py/quickstart.html#primer-design
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER'   : True,
    #'PRIMER_PICK_INTERNAL_OLIGO': True,
    'PRIMER_PICK_RIGHT_PRIMER'  : True,
    'PRIMER_PRODUCT_SIZE_RANGE': [700, 1000],
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 36,
    'PRIMER_OPT_TM': 62.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 88.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 10,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
  }
  # log(args, 'primer_args: %s' % primer_args)
  # log(args, 'seq_args: %s' % seq_args)
  primers = primer3.designPrimers(seq_args, primer_args)
  # TODO parse this into a list of simpler dicts and return thoes
  return primers

# TODO does this mean I should switch to Gibson because it can use almost any sequence equally?

def suggest_left_flank_primers(args, seq):
  return suggest_primers(args, seq, {'SEQUENCE_FORCE_RIGHT_START': len(seq)-1}) # TODO is this right?

def suggest_right_flank_primers(args, seq):
  return suggest_primers(args, seq, {'SEQUENCE_FORCE_LEFT_START': 0})

def extract_result(args, p3res, n):
  'extract info for just one pair of primers from primer3 result dict'
  # TODO return a dict here too, just a simpler one
  # TODO and the numbers can stay floats?
  return (          p3res['PRIMER_LEFT_%d_SEQUENCE'  % n],
          '%0.1f' % p3res['PRIMER_LEFT_%d_TM'        % n],
                    p3res['PRIMER_RIGHT_%d_SEQUENCE' % n],
          '%0.1f' % p3res['PRIMER_RIGHT_%d_TM'       % n])

def simplify_results(args, p3res):
  'group primer3 results by the primer pair (result 1, 2, 3...)'
  # TODO split into left and right dicts each?
  n = 1
  results = []
  while True:
    # keys = [k for k in p3res.keys() if '_%d_' % n in k]
    #print(p3res)
    keys = \
      { 'PRIMER_LEFT_%d_SEQUENCE'     % n: 'left_seq'
      , 'PRIMER_LEFT_%d_TM'           % n: 'left_tm'
      , 'PRIMER_RIGHT_%d_SEQUENCE'    % n: 'right_seq'
      , 'PRIMER_RIGHT_%d_TM'          % n: 'right_tm'
      , 'PRIMER_PAIR_%d_PRODUCT_SIZE' % n: 'size'
      }
    res = {}
    for k in keys.keys():
      try:
        res[keys[k]] = p3res[k]
      except KeyError:
        return results
    log(args, res)
    # keys = [re.sub('PRIMER_LEFT_%d_'  % n, 'left_' , k) for k in keys]
    # keys = [re.sub('PRIMER_RIGHT_%d_' % n, 'right_', k) for k in keys]
    # keys = [re.sub('PRIMER_PAIR_%d_'  % n, 'right_', k) for k in keys]
    #print(keys)
    #print()
    results.append(res)
    n += 1


#######################################################
# add restriction sites to flank primers and re-score #
#######################################################

def add_restriction_sites(args, pairs):
  # kpnI = 'XXXXX'
  # salI = 'YYYYY'
  # TODO reverse complement the rev one?
  # pairs = [(kpnI + fwd
  pass


########
# main #
########

def parse(args):
  return {
    'verbose'  : args['-v'],
    'genome'   : args['-g'],
    'locusids' : ['Synpcc7942_%s' % a for a in args['LOCUS']]
  }

def hr_primers(args, genome, seq, homology_bp=1000):
  locus = seq['id']

  # get flanking sequences
  # left_start  = seq['start'] - homology_bp
  # left_end    = seq['start']
  # right_start = seq['end']
  # right_end   = seq['end'] + homology_bp

  # design primers for them
  # TODO deal with the whole list instead
  left_res  = suggest_primers(args, 'left' , str(genome.seq[left_start  :  left_end]))
  right_res = suggest_primers(args, 'right', str(genome.seq[right_start : right_end]))
  left_fwd_seq , left_fwd_tm , left_rev_seq , left_rev_tm  = extract_result(args, left_res , 1)
  right_fwd_seq, right_fwd_tm, right_rev_seq, right_rev_tm = extract_result(args, right_res, 1)

  # TODO show coordinates separately and just primer stuff here
  log(args, '%s left  flank (%s-%sbp) forward primer: %s (Tm=%s)' % (locus, left_start , left_end , left_fwd_seq , left_fwd_tm ), 2)
  log(args, '%s left  flank (%s-%sbp) reverse primer: %s (Tm=%s)' % (locus, left_start , left_end , left_rev_seq , left_rev_tm ), 2)
  log(args, '%s right flank (%s-%sbp) forward primer: %s (Tm=%s)' % (locus, right_start, right_end, right_fwd_seq, right_fwd_tm), 2)
  log(args, '%s right flank (%s-%sbp) reverse primer: %s (Tm=%s)' % (locus, right_start, right_end, right_rev_seq, right_rev_tm), 2)

  return [left_fwd_seq,  left_fwd_tm,  left_rev_seq,  left_rev_tm,
         right_fwd_seq, right_fwd_tm, right_rev_seq, right_rev_tm]

def main():
  args = parse(docopt(__doc__, version='cpf1primers 0.1'))
  log(args, 'args: %s' % args, 2)
  genome = list(SeqIO.parse(args['genome'], 'genbank'))[0] # TODO look up the right way
  log(args, 'genome: %s' % str(genome.seq)[:50] + '...', 2)

  # this does dicts already
  seqs = find_sequences(args, genome, args['locusids'], homology_bp=1000)

  # the rest, not yet
  # headers = ['locus_ID', 'crRNA_fwd', 'crRNA_rev',
  #             'left_fwd_seq',  'left_fwd_tm', ' left_rev_seq',  'left_rev_tm',
  #            'right_fwd_seq', 'right_fwd_tm', 'right_rev_seq', 'right_rev_tm']
  # print('\t'.join(headers))

  # TODO this should convert each locus dict to a dict of primer dicts:
  #        { 'id': ..., 'crRNA': ..., 'left_flank': {'fwd': ..., 'rev': ....}, ...}
  #        then can add whatever info needed in json, like expected product size and Tm
  primers_by_seq = {}
  for seq in seqs:
      log(args, 'seq: %s' % shorten_seqs(args, seq), 1)
      # row = [seq['id']]
      primers = {}
      primers['crRNA'] = crRNA(args, seq)
      # primers[ 'left_flank'] =  suggest_left_flank_primers(args, seq['left' ])
      # primers[ 'left_flank'] =  suggest_primers(args, seq['left' ])
      primers['left_flank' ] = simplify_results(args,  suggest_left_flank_primers(args, seq['left' ]))
      primers['right_flank'] = simplify_results(args, suggest_right_flank_primers(args, seq['right']))
      # row += [guide_fwd, guide_rev]
      # left_fwd, left_rev, right_fwd, right_rev = hr_primers(args, genome, seq)
      # row += hr_primers(args, genome, seq)
      # print('\t'.join(row))
      log(args, 'primers: %s' % primers, 2)
      primers_by_seq[seq['id']] = primers
  log(args, 'primers by seq: %s' % primers_by_seq, 1)

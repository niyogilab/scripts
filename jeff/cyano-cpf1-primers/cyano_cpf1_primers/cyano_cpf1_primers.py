#!/usr/bin/env python3

# TODO redesign to more generically join pieces 1, 2, 3 or 1, 2, 3, 4 by gibson (are there specific designers for this?)
# TODO 3) add Tms, lengths to output
# TODO 4) redesign for just KpnI digest of backbone, then 3-part gibson
# TODO 5) add some common-sense checks so you don't have to worry about it
# TODO 6) check pre-crRNA design against paper!
# TODO 7) force design of simple gibson primers when primer3 won't cooperate and see if they work
# TODO 8) add NSI knock-ins too so you can order both sets today

# TODO hey you can primer3.calcTm() to check that the overlaps have reasonable Tm on their own!
# TODO add GC_CLAMP?

'''
Generates a table of initial primers to try for Cpf1 knockouts + knock-ins of cyano genes.
Designed for markerless KOs with pSL2680 from the Pakrasi lab (AddGene plasmid #85581).
It defaults to the PCC 7942 genome, but theoretically should work in at least 2973, 7942, 6803, and 7201.
The pre-crRNA pair should be annealed and ligated into the plasmid directly;
left and right flank primers are for generating a standard HR template by PCR with the target gene missing.
The final HR template can be integrated into pSL2680 at the KpnI site.
Complementation primers amplify the target coding sequence and add flanks for the NSI plasmid (???).

Usage:
  cyano_cpf1_primers (-h | --help)
  cyano_cpf1_primers [-v...] [-g GENOME] [-f FORMAT] [-n NOPTS] LOCUS...

Options:
  -h, --help  Show this text
  -v          Print debugging information to stderr
  -g GENOME   PCC 7942 genome to use. [default: SynPCC7942_chr.gbk]
  -f FORMAT   How to format the output. Current options are 'json' or 'table'. [default: json]
  -n NOPTS    Design up to this many alternatives for each primer. [default: 1]
  LOCUS       Locus ID to generate primers for. You can put more than one, and use shorthand.
              For example, 0001 1021 means Synpcc7942_0001 and Synpcc7942_1021.
'''

from __future__  import print_function

import json
import primer3
import re

from Bio          import SeqIO
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from copy         import deepcopy
from docopt       import docopt
from sys          import stderr

# TODO any reason to make this configurable? maybe later
KPNI_SITE = 'GGTACC'

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

def check_no_kpnI_site(args, seq):
  kpnI = KPNI_SITE.lower()
  return len(re.findall(kpnI, seq.lower())) == 0


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

def pre_crRNA_primers_for_cpf1_target(target):
  'generates fwd and rev primers with AarI sites added'
  return {
    'fwd': 'AGAT' + target,
    'rev': 'AGAC' + str(Seq(target, generic_dna).reverse_complement())
  }

def pre_crRNA(args, seq):
  locus = seq['id']
  # sequence = get_sequence(args, locus)
  targets = find_cpf1_targets(args, locus, seq['cds'])
  primers = [pre_crRNA_primers_for_cpf1_target(t) for t in targets[:args['nopts']]]

  # TODO any reason to be more selective?
  # TODO option to pick more than one per locus?
  # TODO can i use primer3 here to to guess which ones will be more effective?
  # target  = targets[0]
  log(args, 'using the first %s: %s' % (args['nopts'], primers), 2)
  # fwd_primers = ['AGAT' + t for t in targets]
  # rev_primers = ['AGAC' + str(Seq(t, generic_dna).reverse_complement()) for t in targets]
  # pairs = zip(fwd_primers, rev_primers)

  # this finds the lowest-Tm pair, but I'm not sure if being nearest the flank is more important
  #(fwd_primer, rev_primer, tm) = pick_lowest_tm(args, pairs))
  # log(args, '%s guide rna annealing temp: %s' % (locus, tm))
  # (fwd_primer, rev_primer) = pairs[0]
  # fwd_primer = primers['fwd']
  # rev_primer = primers['rev']

  log(args, '%s pre-crRNA primers: %s' % (locus, primers), 2)
  return primers


######################################
# generate primer pairs with primer3 #
######################################

# TODO fiddle with these or find sane defaults on the web
# TODO should i have a hybridization probe?
# TODO extract and return just sequences first
# TODO extract Tm too, and add to the table
def suggest_primers(args, seq, extra_seq_args={}, extra_primer_args={}):
  seq_args = {
    'SEQUENCE_ID': 'seqid', # TODO any reason to have use this?
    'SEQUENCE_TEMPLATE': seq,
  }
  seq_args.update(extra_seq_args)
  primer_args = {
    # see https://libnano.github.io/primer3-py/quickstart.html#primer-design
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER'   : True,
    'PRIMER_PICK_INTERNAL_OLIGO': False, # TODO do I want this too for checking?
    'PRIMER_PICK_RIGHT_PRIMER'  : True,
    'PRIMER_PRODUCT_SIZE_RANGE': [700, 1000],

    # TODO how did i decide these?
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
  primer_args.update(extra_primer_args)
  # log(args, 'primer_args: %s' % primer_args)
  # log(args, 'seq_args: %s' % seq_args)
  primers = primer3.designPrimers(seq_args, primer_args)
  # TODO parse this into a list of simpler dicts and return thoes
  return primers

def check_primers(args, seq, left_primer, right_primer):
  seq_args = {
    'SEQUENCE_ID': 'seqid', # TODO any reason to have use this?
    'SEQUENCE_TEMPLATE': seq,
    'SEQUENCE_PRIMER': left_primer,
    'SEQUENCE_PRIMER_REVCOMP': right_primer, # TODO reverse complement it?
  }
  # seq_args.update(extra_seq_args)
  primer_args = {
    # see https://libnano.github.io/primer3-py/quickstart.html#primer-design
    'PRIMER_TASK': 'check_primers',
    'PRIMER_PICK_LEFT_PRIMER'   : False,
    'PRIMER_PICK_INTERNAL_OLIGO': False,
    'PRIMER_PICK_RIGHT_PRIMER'  : False,
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
  primers = primer3.designPrimers(seq_args, primer_args)
  return primers

def assert_primer3(args, seqid, seq, piece, seq_args, primer_args, res):
    if res['PRIMER_LEFT_NUM_RETURNED'] == 0 or res['PRIMER_RIGHT_NUM_RETURNED'] == 0:
        err = '\n'.join([
          'Primer3 found no primers for %s %s' % (seqid, piece),
          'sequence: %s' % seq,
          'seq_args: %s' % seq_args,
          'primer_args: %s' % primer_args,
          'results: %s' % res
        ])
        raise SystemExit(err)

def suggest_left_flank_primers(args, seq, max_th=47, no_thermo=False):
  seq_args = {'SEQUENCE_FORCE_RIGHT_START': len(seq['left'])-1}
  primer_args = {'PRIMER_MIN_SIZE': 2, 'PRIMER_MAX_SIZE': 36, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 0}
  if no_thermo:
    primer_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
  res = suggest_primers(args, seq['left'], seq_args, primer_args) # TODO is this right?
  try:
    assert_primer3(args, seq['id'], seq['left'], 'left', seq_args, primer_args, res)
  except:
    # almost all the failures are due to 'high hairpin stability', so
    # i try relaxing that until it works anyway
    if max_th > 80: # TODO what's a reasonable number?
      if no_thermo:
        raise
      else:
        primer_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
        return suggest_left_flank_primers(args, seq, no_thermo=True)
    else:
      return suggest_left_flank_primers(args, seq, max_th=max_th+1)
  return res

def suggest_right_flank_primers(args, seq, max_th=47, no_thermo=False):
  seq_args = {'SEQUENCE_FORCE_LEFT_START': 0}
  primer_args = {'PRIMER_MIN_SIZE': 2, 'PRIMER_MAX_SIZE': 36, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 0}
  if no_thermo:
    primer_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
  res = suggest_primers(args, seq['right'], seq_args, primer_args)
  # try:
  try:
    assert_primer3(args, seq['id'], seq['right'], 'right', seq_args, primer_args, res)
  except:
    # almost all the failures are due to 'high hairpin stability', so
    # i try relaxing that until it works anyway
    if max_th > 80: # TODO what's a reasonable number?
      if no_thermo:
        rais
      else:
        primer_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
        return suggest_right_flank_primers(args, seq, no_thermo=True)
    else:
      return suggest_right_flank_primers(args, seq, max_th=max_th+1)
  return res

def suggest_coding_primers(args, seq, max_th=47, no_thermo=False):
  seq_args = {
    'SEQUENCE_FORCE_LEFT_START': 0,
    'SEQUENCE_FORCE_RIGHT_START': len(seq['cds'])-1,
  }
  primer_args = {
    'PRIMER_PRODUCT_SIZE_RANGE': [36, 99999], # both ends fixed anyway
    'PRIMER_MIN_SIZE': 2,  # must be > min product size
    'PRIMER_MAX_SIZE': 36, # the built-in max
    'PRIMER_MAX_HAIRPIN_TH': max_th
  }
  if no_thermo:
    primer_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
  res = suggest_primers(args, seq['cds'], seq_args, primer_args)
  try:
    assert_primer3(args, seq['id'], seq['cds'], 'cds', seq_args, primer_args, res)
  except:
    # almost all the failures are due to 'high hairpin stability', so
    # i try relaxing that until it works anyway
    if max_th > 80: # TODO what's a reasonable number?
      if no_thermo:
        raise
      else:
        return suggest_coding_primers(args, seq, no_thermo=True)
    else:
      return suggest_coding_primers(args, seq, max_th=max_th+1)
  return res

def extract_result(args, p3res, n):
  'extract info for just one pair of primers from primer3 result dict'
  # TODO return a dict here too, just a simpler one
  # TODO and the numbers can stay floats?
  return (          p3res['PRIMER_LEFT_%d_SEQUENCE'  % n],
          '%0.1f' % p3res['PRIMER_LEFT_%d_TM'        % n],
                    p3res['PRIMER_RIGHT_%d_SEQUENCE' % n],
          '%0.1f' % p3res['PRIMER_RIGHT_%d_TM'       % n])

def simplify_results(args, p3res):
  'group primer3 results by the primer pair (result 1, 2, 3...) and keep only relevant keys'
  # TODO split into left and right dicts each?
  n = 1
  results = []
  while True:
    # keys = [k for k in p3res.keys() if '_%d_' % n in k]
    #print(p3res)
    keys = \
      { 'PRIMER_LEFT_%d_SEQUENCE'     % n: 'left_primer'
      , 'PRIMER_LEFT_%d_TM'           % n: 'left_tm'
      , 'PRIMER_RIGHT_%d_SEQUENCE'    % n: 'right_primer'
      , 'PRIMER_RIGHT_%d_TM'          % n: 'right_tm'
      , 'PRIMER_PAIR_%d_PRODUCT_SIZE' % n: 'product_size'
      }
    res = {}
    for k in keys.keys():
      try:
        res[keys[k]] = p3res[k]
      except KeyError:
        return results
    log(args, 'primer3 result: %s' % p3res)
    log(args, res)
    # keys = [re.sub('PRIMER_LEFT_%d_'  % n, 'left_' , k) for k in keys]
    # keys = [re.sub('PRIMER_RIGHT_%d_' % n, 'right_', k) for k in keys]
    # keys = [re.sub('PRIMER_PAIR_%d_'  % n, 'right_', k) for k in keys]
    #print(keys)
    #print()
    results.append(res)
    n += 1



########
# main #
########

def parse(args):
  return {
    'verbose'  : args['-v'],
    'genome'   : args['-g'],
    'format'   : args['-f'],
    'nopts'    : int(args['-n']),
    'locusids' : ['Synpcc7942_%s' % a for a in args['LOCUS']] #TODO generalize this
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

def print_row(*vals):
    print('\t'.join(vals))

def print_tsv(args, primers):
    # this is done the stupid way. the smart way would probably involve pandas?
    print('\t'.join(['locusid', 'construct', 'piece', 'primer', 'pair', 'sequence']))
    for locusid in sorted(primers.keys()):
        n = 1
        for pair in primers[locusid]['knockout']['pre_crRNA']:
            print_row(locusid, 'knockout', 'pre-crRNA', 'fwd', str(n), pair['fwd'])
            print_row(locusid, 'knockout', 'pre-crRNA', 'rev', str(n), pair['rev'])
            n += 1
        n = 1
        for pair in primers[locusid]['knockout']['left_flank']:
            print_row(locusid, 'knockout', 'left flank', 'fwd', str(n), pair['left_primer'])
            print_row(locusid, 'knockout', 'left flank', 'rev', str(n), pair['right_primer'])
            n += 1
        n = 1
        for pair in primers[locusid]['knockout']['right_flank']:
            print_row(locusid, 'knockout', 'right flank', 'fwd', str(n), pair['left_primer'])
            print_row(locusid, 'knockout', 'right flank', 'rev', str(n), pair['right_primer'])
            n += 1
        n = 1
        for pair in primers[locusid]['complement']['coding_sequence']:
            print_row(locusid, 'complement', 'coding sequence', 'fwd', str(n), pair[ 'left_primer'])
            print_row(locusid, 'complement', 'coding sequence', 'rev', str(n), pair['right_primer'])
            n += 1

def main():
  args = parse(docopt(__doc__, version='cpf1primers 0.2'))
  log(args, 'args: %s' % args, 2)
  genome = list(SeqIO.parse(args['genome'], 'genbank'))[0] # TODO look up the right way
  log(args, 'genome: %s' % str(genome.seq)[:50] + '...', 2)

  # this does dicts already
  seqs = find_sequences(args, genome, args['locusids'], homology_bp=1000)

  # the rest, not yet
  # headers = ['locus_ID', 'pre_crRNA_fwd', 'pre_crRNA_rev',
  #             'left_fwd_seq',  'left_fwd_tm', ' left_rev_seq',  'left_rev_tm',
  #            'right_fwd_seq', 'right_fwd_tm', 'right_rev_seq', 'right_rev_tm']
  # print('\t'.join(headers))

  primers_by_seq = {}
  for seq in seqs:
      log(args, 'seq: %s' % shorten_seqs(args, seq), 1)
      # row = [seq['id']]

      # generate knockout primers
      ko = {}
      ko['pre_crRNA'] = pre_crRNA(args, seq)
      assert check_no_kpnI_site(args, seq['left'])
      assert check_no_kpnI_site(args, seq['right'])
      log(args, 'cool, no existing kpnI sites in the flanks.')
      # here are the initial primers designed by primer3
      ko['left_flank' ] = simplify_results(args,  suggest_left_flank_primers(args, seq))[:args['nopts']]
      ko['right_flank'] = simplify_results(args, suggest_right_flank_primers(args, seq))[:args['nopts']]
      # TODO now we add restriction sites to each pair
      # TODO along with at least 2 more bp overhang
      # TODO and have primer3 re-score them before deciding on the best? or ipcress or something else
      # TODO WHEN ADDING GIBSON OVERLAPS, NEED TO DEAL WITH PAIRS OF PRIMER PAIRS INSTEAD OF EACH SEPARATE
      pairs2 = []
      for pair in ko['left_flank']:
        l = KPNI_SITE.lower() + pair['left_primer']
        r = str(Seq(l[-15:], generic_dna).reverse_complement()).lower() + pair['right_primer']
        pair.update({'left_primer': l, 'right_primer': r})
        pairs2.append(pair)
        # pairs2.append((add_restriction_sites(args, seq['left'], pair)))
      ko['left_flank'] = pairs2
      pairs2 = []
      for pair in ko['right_flank']:
        r = KPNI_SITE.lower() + pair['right_primer']
        l = str(Seq(r[-15:], generic_dna).reverse_complement()).lower() + pair['left_primer']
        pair.update({'left_primer': l, 'right_primer': r})
        pairs2.append(pair)
      ko['right_flank'] = pairs2
      # row += [guide_fwd, guide_rev]
      # left_fwd, left_rev, right_fwd, right_rev = hr_primers(args, genome, seq)
      # row += hr_primers(args, genome, seq)
      # print('\t'.join(row))
      log(args, 'ko: %s' % ko, 2)

      # TODO generate knock-in primers
      cds = simplify_results(args, suggest_coding_primers(args, seq))[:args['nopts']]
      # print(cds)
      ki = {'coding_sequence': cds}
      log(args, 'ki: %s' % ki, 2)
      # TODO always use the same NSI pre-crRNA here? or skip if the pakrasi paper provides it
      # TODO suggest left/right for coding sequence (will this usually fail??)
      # TODO add overlap from the plasmid backbone (hardcoded first, then configurable later)

      primers_by_seq[seq['id']] = {'knockout': ko, 'complement': ki}

  if args['format'] == 'table':
    print_tsv(args, primers_by_seq)
  else:
    log(args, json.dumps(primers_by_seq, sort_keys=True, indent=2), 0)

from Bio import SeqIO
from collections import defaultdict
from itertools import permutations
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, nargs='+',
	help='fasta files to be concatenated')
parser.add_argument('-d','--hmmdir', type=str,
	help='Directory containing hmm output')
args = parser.parse_args()

def parse_fasta(infile):
    return {rec.id:rec for rec in SeqIO.parse(infile, 'fasta')}

def get_domains(infile):
    domains = defaultdict(dict)
    count = defaultdict(int)
    for line in open(infile):
        if not line.startswith('#'):
            line = line.strip().split()
            recid = line[0]
            domain = line[3]
            tstart = line[17]
            tend = line[18]
            score = line[7]
            if domain in domains[recid]:
                count[recid] += 1
                domain = domain + '_{}'.format(count[recid])
            domains[recid][domain] = {'tstart':int(tstart), 'tend':int(tend), 'score':float(score), 'length':int(tend)-int(tstart)}
    return domains

def resolve_overlaps(d, b):
    if (d['tend'] >= b['tstart'] and d['tend'] <= b['tend']):
        if d['score'] > b['score']:
            b['tstart'] = d['tend'] + 1
        elif b['score'] > d['score']:
            d['tend'] = b['tstart'] - 1
    elif (d['tstart'] <= b['tend'] and d['tstart'] >= b['tstart']):
        if d['score'] > b['score']:
            b['tend'] = d['tstart'] - 1
        elif b['score'] > d['score']:
            d['tstart'] = b['tend'] + 1

def resolve_incompatible_domains(domains):
    for k,v in domains.items():
        for (a,b),(c,d) in permutations(v.items(),2):
            resolve_overlaps(b,d)
    for k,v in domains.items():
        for a,b in list(v.items()):
            if b['tend'] - b['tstart'] <= .3 * b['length']:
                v.pop(a)

def exctract_subseq(seq, domain, meta):
    subseq = seq[meta['tstart']:meta['tend']]
    subseq.id = seq.id + '_' + domain
    subseq.description = ''
    return subseq

def main(base, hmmdir):
    seqdict = parse_fasta(fasta)
    base = os.path.basename(fasta).replace('.fasta', '')
    domains = get_domains('{}/{}.out'.format(hmmdir,base))
    resolve_incompatible_domains(domains)

    for split_file in set([k.split('_')[0] for i in domains.values() for k in i.keys()]):
        subseqs = []
        for seq, doms in domains.items():
            for dom, meta in doms.items():
                if dom.split('_')[0] == split_file:
                    subseqs.append(exctract_subseq(seqdict[seq], dom, meta))
        if len(subseqs) >= 0.0 * len(seqdict.values()):
            with open("COG_split_resolve_eq/{}_{}.fasta".format(base, split_file), 'w') as out:
                SeqIO.write(subseqs, out, 'fasta')

if __name__ == '__main__':
    for fasta in args.fasta:
        main(fasta, args.hmmdir)

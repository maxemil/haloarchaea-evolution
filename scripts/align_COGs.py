import glob
import os
from multiprocessing import Pool
import subprocess
from collections import defaultdict
from Bio import SeqIO


eur2NOG = {line.split()[0]:line.strip().split()[1].split(';') for line in open('inconsistent_eurNOGs.txt')}
cogs = list(set([cog for l in eur2NOG.values() for cog in l]))
#
# cog2prot = defaultdict(list)
# for line in open('cog2003-2014.csv'):
#     line = line.strip().split(',')
#     if line[6] in cogs:
#         cog2prot[line[6]].append([line[0], line[4], line[5]])
#
# seqdict = {rec.id.split('|')[1]: rec for rec in SeqIO.parse('prot2003-2014.fa', 'fasta')}
#
# for cog in cogs:
#     with open('COG/{}.fasta'.format(cog), 'w') as out:
#         for prot in cog2prot[cog]:
#             rec = seqdict[prot[0]][int(prot[1])-1:int(prot[2])]
#             SeqIO.write(rec, out, 'fasta')
#
# def cluster_sequences(fasta):
#     subprocess.check_call('cd-hit -i {} -o {} -c 0.8'.format(fasta, fasta + ".reduced").split())
#
# p = Pool(30)
# p.map(cluster_sequences, glob.glob('COG/*.fasta'))
#
# def align_hmm(fasta):
#     if not (os.path.isfile(fasta.replace('.fasta.reduced', '.aln'))):
#         subprocess.check_call('mafft-fftnsi --anysymbol --thread 6 {} '.format(fasta).split(), stdout=open(fasta.replace('.fasta.reduced', '.aln'), 'w'))
#     if not (os.path.isfile(fasta.replace('.fasta.reduced', '.hmm'))):
#         subprocess.check_call('hmmbuild {} {}'.format(fasta.replace('.fasta.reduced', '.hmm'), fasta.replace('.fasta.reduced', '.aln')).split())
#
# p = Pool(7)
# p.map(align_hmm, glob.glob('COG/*.fasta.reduced'))

def concat_hmmsearch(fasta):
    base = fasta.replace('clusters_to_split/', '').replace('.fasta', '')
    refs = ['COG/{}.hmm'.format(nog) for nog in eur2NOG[base] if os.path.isfile('COG/{}.hmm'.format(nog))]
    subprocess.check_call(["cat"] + refs, stdout=open("HMM_output/{}.hmm".format(base), 'w'))
    subprocess.check_call("hmmsearch -o {} --max --domtblout {} --domE 0.001 --cpu 1 {} {}".format("HMM_output/{}.log".format(base), "HMM_output/{}.out".format(base), "HMM_output/{}.hmm".format(base), fasta).split())

p = Pool(2)
p.map(concat_hmmsearch, glob.glob('clusters_to_split/*'))

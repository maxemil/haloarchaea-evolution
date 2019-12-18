import glob
import os
from multiprocessing import Pool
import subprocess


# cog2prot = defaultdict(list)
# for line in open('COG/cog2003-2014.csv'):
#     line = line.strip().split(',')
#     if line[6] in cogs:
#         cog2prot[line[6]].append([line[0], line[4], line[5]])
#
# seqdict = {rec.id.split('|')[1]: rec for rec in SeqIO.parse('COG/prot2003-2014.fa', 'fasta')}
#
# for cog in cogs:
#     with open('COG/{}.fasta'.format(cog), 'w') as out:
#         for prot in cog2prot[cog]:
#             rec = seqdict[prot[0]][int(prot[1])-1:int(prot[2])]
#             SeqIO.write(rec, out, 'fasta')


def align_hmm(fasta):
    # if not (os.path.isfile(fasta.replace('.fasta.reduced', '.aln'))):
        # subprocess.check_call('mafft-fftnsi --anysymbol --thread 10 {} '.format(fasta).split(), stdout=open(fasta.replace('.fasta.reduced', '.aln'), 'w'))
    subprocess.check_call('hmmbuild {} {}'.format(fasta.replace('.fasta.reduced', 'hmm'), fasta.replace('.fasta.reduced', '.aln')).split())

def concat_hmmsearch(fasta):
    base = fasta.replace('clusters/', '').replace('.fasta', '')
    # refs = ['COG/{}.hmm'.format(nog) for nog in eur2NOG[base] if os.path.isfile('COG/{}.hmm'.format(nog))]
    # subprocess.check_call(["cat"] + refs, stdout=open("COG/{}.hmm".format(base), 'w'))
    subprocess.check_call("hmmsearch -o {} --max --domtblout {} --domE 0.001 --cpu 5 {} {}".format("COG/{}.log".format(base), "COG/{}.out".format(base), "COG/{}.hmm".format(base), fasta).split())

# p = Pool(8)
# p.map(align_hmm, glob.glob('COG/*.fasta.reduced'))


eur2NOG = {line.split()[0]:line.strip().split()[1].split(';') for line in open('../../20_clusters_eggnog/inconsistent_eurNOGs.txt')}
p = Pool(8)
p.map(concat_hmmsearch, glob.glob('clusters/*'))

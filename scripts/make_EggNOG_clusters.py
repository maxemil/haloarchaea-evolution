import glob
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import argparse
import requests
import os
import sys
import shutil
from subprocess import Popen

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', type=str,
            help='glob to fasta files containing protein sequences')
parser.add_argument('-a', '--annotations', type=str,
            help='glob to annotation files from emapper, matching fasta file sequences')
parser.add_argument('-o', '--outdir', type=str,
            help='output_directory')
parser.add_argument('-s', '--hifix', type=str,
            help='location of HiFiX.simg')
parser.add_argument('-r', '--references', action='store_true', help="download reference sequences")
args = parser.parse_args()

def parse_annotations(annotations):
    eggnog_dict = defaultdict(list)
    for annot in glob.glob(annotations):
        df = pd.read_csv(annot, skiprows=4, skipfooter=3, sep='\t', index_col=False, header=None, engine='python')
        df = df.apply(lambda row:  eggnog_dict[row[10].split('|')[0]].append(row[0]), axis=1)
    return eggnog_dict

def parse_fasta(fasta):
    seq_dict = {}
    for fasta in glob.glob(fasta):
        for rec in SeqIO.parse(fasta, 'fasta'):
            seq_dict[rec.id] = rec
    return seq_dict

def write_clusters(eggnog_dict, seq_dict, outdir):
    for nog in eggnog_dict.keys():
        with open("{}/{}.fasta".format(outdir, nog), 'w') as out:
            for seq in eggnog_dict[nog]:
                SeqIO.write(seq_dict.pop(seq), out, 'fasta')
    return seq_dict

def write_unlabelled_seqs(seq_dict, outdir):
    with open('{}/silix/seqs.fasta'.format(outdir), 'w') as out:
        for seq in list(seq_dict.keys()):
            SeqIO.write(seq_dict.pop(seq), out, 'fasta')
    return seq_dict

def cluster_unlabelled_seqs(outdir, hifix):
    p = Popen("diamond makedb --in {}/silix/seqs.fasta --db {}/silix/seqs.dmnd".format(outdir, outdir).split())
    p.communicate()
    p = Popen("diamond blastp --query {}/silix/seqs.fasta --db {}/silix/seqs.dmnd --out {}/silix/seqs.blastp --threads 40 --outfmt 6 --max-target-seqs 1000 --more-sensitive --evalue 0.0001".format(outdir, outdir, outdir).split())
    p.communicate()
    with open("{}/silix/seqs.blastpath".format(outdir), 'w') as out:
        p = Popen("echo {}/silix/seqs.blastp".format(outdir).split(), stdout=out)
        p.communicate()
    with open("{}/silix/unlab_60_90.fnodes".format(outdir), 'w') as out:
        p = Popen("singularity exec -B {} {} silix {}/silix/seqs.fasta {}/silix/seqs.blastpath --prefix 60_90_ --ident 0.6 --overlap 0.9 --net".format(outdir, hifix, outdir, outdir).split(), stdout=out)
        p.communicate()
    p = Popen("singularity exec -B {} {} silix-split -n 1 {}/silix/seqs.fasta {}/silix/unlab_60_90.fnodes -p {}/unlab".format(outdir, hifix, outdir, outdir, outdir).split())
    p.communicate()


def download_from_eggnog(attribute, nogname, gzip, outdir):
    url = "http://eggnogapi.embl.de/nog_data/file/{}/{}".format(attribute, nogname)
    r = requests.get(url, stream=True)
    if gzip:
        filename = "{}/{}.{}.gz".format(outdir, nogname, attribute)
    else:
        filename = "{}/{}.{}".format(outdir, nogname, attribute)
    with open(filename, 'wb') as f:
        shutil.copyfileobj(r.raw, f)

def main(args):
    if os.path.exists(args.outdir):
        sys.exit("outdir {} already exists, exiting...".format(args.outdir))
    else:
        os.makedirs(args.outdir)
        os.makedirs("{}/silix".format(args.outdir))
    eggnog_dict = parse_annotations(args.annotations)
    seq_dict = parse_fasta(args.fasta)
    seq_dict = write_clusters(eggnog_dict, seq_dict, args.outdir)
    seq_dict = write_unlabelled_seqs(seq_dict, args.outdir)
    assert not seq_dict
    cluster_unlabelled_seqs(args.outdir, args.hifix)
    if args.references and not os.path.exists("references"):
        os.makedirs("references")
        for eggnog in eggnog_dict.keys():
            download_from_eggnog('fasta', eggnog, True, "references")
            download_from_eggnog('tree', eggnog, False, "references")
            download_from_eggnog('raw_alg', eggnog, True, "references")

if __name__ == '__main__':
    main(args)

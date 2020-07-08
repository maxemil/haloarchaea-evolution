from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', type=str,
            help='input fasta file')
parser.add_argument('-o', '--output', type=str,
            help='output fasta file')
parser.add_argument('-p', '--prefix', type=str,
            help='prefix to add to sequences')
parser.add_argument('-s', '--separator', type=str, default='..',
            help='separator between species and sequence name')
args = parser.parse_args()

with open(args.output, 'w') as outfile, open(args.fasta + '.log' ,'w') as logfile:
    for rec in SeqIO.parse(args.fasta, 'fasta'):
        new_id = "{}{}{}".format(args.prefix, args.separator, rec.id)
        print("{}\t{}".format(rec.id, new_id), file=logfile)
        rec.id = new_id
        rec.description = ""
        SeqIO.write(rec, outfile, 'fasta')

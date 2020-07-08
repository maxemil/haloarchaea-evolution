"""find inconsistent NOGs

Usage:
  find_inconsistent_NOGs.py -a ANNOTATIONS
  find_inconsistent_NOGs.py --annnotations ANNOTATIONS
  find_inconsistent_NOGs.py (-h | --help)

Options:
  -h --help                     Show this screen.
  -a --annotations ANNOTATIONS  Directory containing .annotations files from Eggnog-mapper
"""
from collections import defaultdict
import glob
import pandas as pd
from docopt import docopt

def get_cog_nog(row, eggnog_dict):
    for l in row[9].split(','):
        l = l.split('@')
        if l[1] == 'NOG':
            eggnog_dict[row[10].split('|')[0]][l[0]] += 1

def int_defaultdict():
    return defaultdict(int)

def main(annotations):
    eggnog_dict = defaultdict(int_defaultdict)
    for annot in glob.glob("{}/*.annotations".format(annotations)):
        df = pd.read_csv(annot, skiprows=4, skipfooter=3, sep='\t', index_col=False, header=None, engine='python')
        df.apply(lambda row: get_cog_nog(row, eggnog_dict) , axis=1)

    for k, v in eggnog_dict.items():
        if len(v) > 1:
            print(k + "\t" +  ";".join(list(v)))

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args['--annotations'])

from collections import defaultdict
import glob
import pandas as pd

def get_cog_nog(row, eggnog_dict):
    for l in row[9].split(','):
        l = l.split('@')
        if l[1] == 'NOG':
            eggnog_dict[row[10].split('|')[0]][l[0]] += 1

def int_defaultdict():
    return defaultdict(int)

eggnog_dict = defaultdict(int_defaultdict)
for annot in glob.glob("../annotations_4.5.1_eurNOG/*.annotations"):
    df = pd.read_csv(annot, skiprows=4, skipfooter=3, sep='\t', index_col=False, header=None, engine='python')
    df.apply(lambda row: get_cog_nog(row, eggnog_dict) , axis=1)

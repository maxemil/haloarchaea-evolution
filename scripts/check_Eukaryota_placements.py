import pandas as pd
from scipy import stats
from ete3 import ncbi_taxonomy
ncbi = ncbi_taxonomy.NCBITaxa()
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO
from statannot import add_stat_annotation

def check_lineage(taxid):
    haleuk = ['Halobacteria', 'Eukaryota']
    lineage = ncbi.get_lineage(int(taxid))
    try:
        l = list(ncbi.get_taxid_translator(lineage).values())
        if any([he in l for he in haleuk]):
            return True
        else:
            return False
    except:
        print("failed to open {}".format(taxid))
        return False

df = pd.read_csv('all_placements.csv', sep='\t')
df = df[df['COG'] != 'unlab']
df_euk = df[df['group'] == "Eukaryota"]

cog2length = {}
cog2tax = {}
for cog in df_euk['COG']:
    if not cog in cog2length.keys():
        try:
            aln = AlignIO.read("references_trimal/{}.trimmed_alg".format(cog), 'fasta')
            cog2length[cog] = (len(aln), aln.get_alignment_length())
            if all([check_lineage(rec.id.split('.')[0]) for rec in aln]):
                cog2tax[cog] = "HalEuk"
            else:
                cog2tax[cog] = "diverse"
        except:
            print("failed to open {}".format(cog))

for cog in df['COG']:
    if not cog in cog2length.keys():
        try:
            aln = AlignIO.read("references_trimal/{}.trimmed_alg".format(cog), 'fasta')
            cog2length[cog] = (len(aln), aln.get_alignment_length())
        except:
            print("failed to open {}".format(cog))

euks_tax = []
euks_len = []
rest_tax = []
rest_len = []
for line in df.iterrows():
    cog = line[1]['COG']
    if cog in cog2length.keys():
        if line[1]['group'] == 'Eukaryota':
            euks_tax.append(cog2length[cog][0])
            euks_len.append(cog2length[cog][1])
        else:
            rest_tax.append(cog2length[cog][0])
            rest_len.append(cog2length[cog][1])

print(stats.ttest_ind(euks_tax, rest_tax, equal_var=False))
print(stats.ttest_ind(euks_len, rest_len, equal_var=False))


def get_length(cog):
    try:
        return cog2length[cog][0]
    except:
        return np.nan

df['length'] = df['COG'].apply(get_length)
plot_order = sorted(set(df['group']), key=lambda c: df[df['group'] == c]['length'].mean(), reverse=True)

plt.figure(figsize=(17,10))
ax = sns.boxplot(x="group", y="length", data=df, order=plot_order)
ax.set(ylabel='Number of proteins in reference NOG')
ax.set(xlabel='Group')
test_results = add_stat_annotation(ax, x="group", y="length", data=df, order=plot_order,
                                   box_pairs=[("Eukaryota", "Proteobacteria"), ("Eukaryota", "FCB group"), ("Eukaryota", "cellular organisms"), ("Eukaryota", "Bacteria"), ("Eukaryota", "others")],
                                   test='t-test_ind', text_format='star',
                                   loc='outside', verbose=2)
plt.tight_layout()
plt.savefig("prot_count_refOG_new.pdf")
plt.close()

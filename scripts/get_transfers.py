from collections import defaultdict
import glob
import re
from math import modf
import seaborn as sns
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import ete3

def get_nodes(tree, clade):
    subtree = tree & clade
    clade_nodes = [n.name for n in subtree.traverse()]
    return clade_nodes

tree = ete3.PhyloTree('species_trees/MethHikHalo_f30_19000_clean.tree')
for n in tree.traverse():
    if not n.is_leaf(): 
        n.name = str(int(n.support))

clades = {
'Halobacteria':'107',
'Hikarchaea':'93',
'Methanomicrobia':'98',
'Archaeoglobales':'101',
'Methanocellales':'91',
'Methanosarcinales':'104',
'Syntrophoarchaea':'62'
}
clades = {v:get_nodes(tree, k) for v,k in clades.items()}
clades = {s: c for c, specs in clades.items() for s in specs}

transfers = defaultdict(lambda: [0, 0])

transfer_nodes = defaultdict(int)

pattern = re.compile(".T@([A-Z0-9]+)->([A-Z0-9]+):")

for umlrec in glob.glob('MGIV-eurNOG/MethHikHalo_f30_19000_clean/*.uml_rec'):
    trans = []
    for line in open(umlrec):
        if line.startswith('('):
            mat = re.findall(pattern, line)
            trans = trans + mat
    for i in set(trans):
        freq = trans.count(i)/100
        freq = modf(freq)[1] + int(modf(freq)[0] >= 0.3)
        transfer_nodes[i] += freq
        if i[1] in set(clades.keys()) and i[0] in set(clades.keys()):
            if clades[i[1]] == clades[i[0]]:
                transfers[i[1]][0] += freq
            else:
                transfers[i[1]][1] += freq

with open('transfers_per_node.csv', 'w') as out:
    print("\t".join(["source", "target", "# transfers > 0.3", "source group", "target group"]), file=out)
    for k, v in transfer_nodes.items():
        if v:
            if k[1] in set(clades.keys()) and k[0] in set(clades.keys()):
                print("\t".join([k[0], k[1], str(v), clades[k[0]], clades[k[1]]]), file=out)
            else:
                print("\t".join([k[0], k[1], str(v), "-", "-"]), file=out)

header = ['species_tree', 'cluster_ID', 'node', 'duplications',
          'transfers', 'losses', 'originations', 'copies']
# df = pd.read_csv("MGIV-eurNOG/events.txt", sep='\t', names=header)
df = pd.read_csv("MGIV-eurNOG/events.txt", sep='\t', names=header)
df['copies'] = df['copies'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
df['transfers'] = df['transfers'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
df['originations'] = df['originations'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
df['losses'] = df['losses'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
df['duplications'] = df['duplications'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))

plt.figure(figsize=(12,8))
freqs = pd.DataFrame()
for k, (i,o) in transfers.items():
    genome_size = df[df['node'] == k]['copies'].sum()
    freqs = freqs.append(pd.DataFrame({"Frequency":o/genome_size, "Clade":clades[k], 'Direction':'intergroup'}, index=[0]))
    freqs = freqs.append(pd.DataFrame({"Frequency":i/genome_size, "Clade":clades[k], 'Direction':'intragroup'}, index=[0]))
    print(k, genome_size, i, o)
ax = sns.swarmplot(x="Clade", y="Frequency", data=freqs, hue='Direction')
ax.set(ylabel='Frequency [normalized with genome size]')
plt.savefig("transfers_new.pdf")
plt.close()

for clade in set(clades.values()):
    clade_freq = freqs[freqs['Clade'] == clade]
    intra = list(clade_freq[clade_freq['Direction'] == 'intragroup']['Frequency'])
    inter = list(clade_freq[clade_freq['Direction'] == 'intergroup']['Frequency'])
    print(clade)
    print(stats.ttest_ind(inter, intra, equal_var=False))

# with open('transfers_03.tab', 'w') as out:
#     for umlrec in glob.glob('MGIV-eurNOG/eurySel.noSA1.chi2prune_clean/*.uml_rec'):
#         trans = []
#         for line in open(umlrec):
#             if line.startswith('('):
#                 mat = re.findall(pattern, line)
#                 trans = trans + mat
#         for i in set(trans):
#             if trans.count(i)/100 >= 0.3:
#                 print(os.path.basename(umlrec).split('.')[0], i[0], i[1], trans.count(i)/100, file=out)

import os
import pandas as pd
import glob
from ete3 import ncbi_taxonomy
from Bio import SeqIO
from collections import defaultdict
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
ncbi = ncbi_taxonomy.NCBITaxa()
sns.set()

REF_TAXA = set(['Proteobacteria', 'Euryarchaeota', 'Terrabacteria group',
                'Acidobacteria', 'Spirochaetes', 'TACK group', 'PVC group',
                'FCB group', 'Deferribacteres', 'Aquificae', 'DPANN group',
                'Synergistetes', 'Eukaryota', 'Thermotogae',
                'Thermodesulfobacteriaceae', 'Spirochaetes', 'Haloplasmatales'])

OTHERS = set(['Acidobacteria', 'Aquificae', 'Deferribacteres',
              'DPANN group', 'PVC group', 'Synergistetes', 'Thermotogae',
              'Thermodesulfobacteriaceae', 'Spirochaetes', 'Haloplasmatales'])

CAT_GROUPS = {"Metabolism":['C','G','E','F','H','I','P','Q', 'C;E', 'E;P'],
              "Cellular Processes":['D','Y','V','T','M','N','Z','W','U','O'],
              "Information":['J','A','K','L','B'],
              "Unknown":['S']}


def get_phylum(taxon):
    if taxon == "Nitrosoarchaeum":
        taxon = "Nitrosarchaeum"
    if taxon == 'Cyanothece sp. ATCC 51142':
        taxon = "Crocosphaera subtropica ATCC 51142"
    try:
        lineage = ncbi.get_lineage(ncbi.get_name_translator([taxon])[taxon][0])
        for l in lineage:
            lname = ncbi.get_taxid_translator([l])[l]
            if lname in REF_TAXA:
                if lname in OTHERS:
                    return 'others'
                return lname
        return taxon
    except:
        return taxon


def parse_placements(placements):
    tax = pd.DataFrame()
    for placement in placements:
        base = os.path.basename(placement).split('.')[0]
        cog = os.path.basename(placement).split('.')[1]
        df = pd.read_csv(placement, sep='\t')
        taxon = df['taxopath'][df['fract'].idxmax()].split(';')[-1]
        tax = tax.append(pd.DataFrame({"COG":cog, "tax":taxon}, index=[base]))
    return tax


def parse_hierarchy(hierarchy_file):
    eur2NOG = {}
    for line in open(hierarchy_file):
        line = line.strip().split()
        nog = line[0]
        for eur in line[1].split(';'):
            eur2NOG[eur] = nog
    return eur2NOG


def get_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(int(taxid))
        lineage_string = ""
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            lineage_string = "{} {};".format(lineage_string, ncbi.get_taxid_translator([l])[l])
        return lineage_string
    except:
        return None


def commonprefix(taxon_paths):
    taxon_path_lists = [t.split('; ') for t in taxon_paths if t is not None]
    s1 = min(taxon_path_lists)
    s2 = max(taxon_path_lists)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1


def get_tax_failed_ref(taxids, ref):
    remaining_taxa = []
    for rec in SeqIO.parse("references_trimal/{}.trimmed_alg".format(ref), 'fasta'):
        taxon = rec.id.split('.')[0]
        if not taxon in taxids:
            remaining_taxa.append(get_lineage(taxon))
    if remaining_taxa:
        return commonprefix(remaining_taxa)[-1].strip(';')
    else:
        return "de-novo genes"


def include_remaining_queries(tax, eur2NOG, taxids):
    for f in glob.glob("queries_all/*"):
        eur = f.replace('queries_all/', '').replace('.fasta', '')
        if eur.startswith('unlab'):
                tax = tax.append(pd.DataFrame({"COG":'Unknown',
                                           "tax":'de-novo genes'},
                                            index=[eur]))
        elif not os.path.exists('placements/{}.{}.csv'.format(eur, eur2NOG[eur])):
            try:
                tax = tax.append(pd.DataFrame({"COG":eur2NOG[eur],
                                           "tax":get_tax_failed_ref(taxids, eur2NOG[eur])},
                                            index=[eur]))
            except:
                print("something else is wrong with reference {} of {}".format(eur2NOG[eur], eur))
            # try:
            #     # if [rec.id.split('.')[0] in taxids for rec in SeqIO.parse("failed_references/{}.raw_alg".format(eur2NOG[eur]),'fasta')].count(False) < 4:
            #     #     tax = tax.append(pd.DataFrame({"COG":eur2NOG[eur], "tax":"de-novo genes"}, index=[eur]))
            # except:
            #     print("something else is wrong with reference {} of {}".format(eur2NOG[eur], eur))
    return tax


def get_clst2node(gain_tables):
    clst2node = defaultdict(lambda: [])
    node2tab = {}
    for t in gain_tables:
        node = t.strip('.tab').split('to')[1]
        df = pd.read_csv(t, sep='\t')
        node2tab[node] = df
        for c in df['cluster']:
            clst2node[c].append(node)
            # if not c in clst2node.keys():
            #     clst2node[c] = node
            # else:
            #     freq_curr = node2tab[node][node2tab[node]['cluster'] == c]['acquisition_freq'].item()
            #     freq_prev = node2tab[clst2node[c]][node2tab[clst2node[c]]['cluster'] == c]['acquisition_freq'].item()
            #     if float(freq_curr) > float(freq_prev):
            #         clst2node[c] = node
            #     else:
            #         pass
            #         # print("{} occurs more than once, currenct is {} at {}, prev was {} at {}".format(c, freq_curr, node, freq_prev, clst2node[c]))
    return(clst2node, node2tab)


def parse_losses(loss_tables):
    loss = pd.DataFrame()
    for t in loss_tables:
        node = t.strip('.tab').split('to')[1]
        n2cat = defaultdict(int)
        for line in open(t):
            if not line.startswith('cluster'):
                cat = line.split('\t')[2]
                n2cat[cat] += 1
        for k,v in n2cat.items():
            loss = loss.append(pd.DataFrame({"cat":k, "freq":v, "node":node}, index=[0]))
    return loss


def get_catmap():
    catmap = defaultdict(lambda: "Unknown")
    for k, value in CAT_GROUPS.items():
        for v in value:
            catmap[v] = k
    return catmap


def print_raw_results(tax, annot, node):
    catmap = get_catmap()
    df = tax.merge(annot, left_index=True, right_index=True, how='left').merge(node, left_index=True, right_index=True, how='left').fillna('Unknown')
    df["event"] = "gain"
    df['group'] = df['tax'].apply(get_phylum)
    df['cat_group'] = df['cat'].apply(lambda x: catmap[x])
    df['major_group'] = df.apply(add_column_major_groups_only, axis=1)
    df.to_csv('all_placements.csv', sep='\t', index=True, index_label="eurNOG" )


def convert_merge_dfs(tax, annot, node, loss):
    catmap = get_catmap()
    annot['cat'] = annot['cat'].apply(lambda x: catmap[x])
    loss['cat'] = loss['cat'].apply(lambda x: catmap[x])
    df = tax.merge(annot, left_index=True, right_index=True, how='left').merge(node, left_index=True, right_index=True, how='left').fillna('Unknown')
    df['tax'] = df['tax'].apply(get_phylum)
    freq = df.groupby(['node', 'tax','cat']).size().reset_index()
    freq.columns = ['node', 'tax', 'cat', 'freq']
    freq["event"] = "gain"
    loss["event"] = "loss"
    freq = pd.concat([freq, loss], sort=True)
    freq = freq[['tax', 'cat', 'event', 'freq', 'node']]
    freq.reset_index(drop=True, inplace=True)
    return freq


def barplot_placements(g, ax, hchars=['//', '.', '|', '\\'], cols=['#332288', '#DDCC77', '#882255', '#44AA99']):
    ordered_tax = ['de-novo genes', 'others', 'cellular organisms', 'Eukaryota',
                    'Archaea', 'TACK group','Euryarchaeota',
                    'Bacteria', 'Proteobacteria', 'Terrabacteria group','FCB group']
    ordered_tax.reverse()
    cmap = ListedColormap(sns.color_palette(cols).as_hex())
    g = g[g['event'] =='gain']
    g = g.drop('node', axis=1)
    g = g.drop('event', axis=1)
    g = g.pivot(index='tax', columns='cat', values='freq')
    g = g.reindex(ordered_tax).fillna(0)#.dropna(how='all')
    g.plot(kind='barh', stacked=True, ax=ax, colormap=cmap, xlim=(0,200))
    # bars = ax.patches
    # hatches = [h for h in hchars for _ in range(len(g))]
    # for bar, hatch in zip(bars, hatches):
    #     bar.set_hatch(hatch)


def plot_ancestors(node_dict):
    fig, axs = plt.subplots(3, figsize=(15,10), sharex=True)
    barplot_placements(node_dict['107'], axs[0])
    barplot_placements(node_dict['93'], axs[1])
    barplot_placements(node_dict['108'], axs[2])
    plt.savefig('ancestors_bars.pdf')
    plt.close()


def plot_all(node_dict):
    fig, axs = plt.subplots(nrows=12, ncols=5, sharex=True, sharey=True, figsize=(30,50))
    for (node, g), ax in zip(node_dict.items(), axs.flatten()):
        # ax = plt.figure(figsize=(20, 10)).add_subplot(111)
        try:
            barplot_placements(g, ax)
            ax.set_title(node)
        except TypeError:
            ax.set_title(node)
    plt.savefig('all_bars_scaleX.pdf')
    plt.close()

others2group = {
'Thermodesulfobacteriaceae':'Bacteria',
'Deferribacteraceae':'Bacteria',
'Aquificae':'Bacteria',
'PVC group':'Bacteria',
'Lentisphaeria':'Bacteria',
'Hydrogenobacter thermophilus TK-6':'Bacteria',
'Nanoarchaeum equitans Kin4-M':'Archaea',
'Pirellula staleyi DSM 6068':'Bacteria',
'Haloplasma contractile SSD-17B':'Bacteria',
'Synergistaceae':'Bacteria',
'Candidatus Solibacter usitatus Ellin6076':'Bacteria',
'Candidatus Koribacter versatilis Ellin345':'Bacteria',
"Leptospira biflexa serovar Patoc strain 'Patoc 1 (Paris)'":'Bacteria',
'Pedosphaera parvula Ellin514':'Bacteria',
'Parachlamydia acanthamoebae UV-7':'Bacteria',
'Planctomycetaceae':'Bacteria',
'Verrucomicrobia':'Bacteria',
'Acidobacteriaceae':'Bacteria',
'Granulicella':'Bacteria',
'Spirochaetia':'Bacteria',
'Thermotogae':'Bacteria'
}

def add_column_major_groups_only(t):
    if t['group'] in ['Proteobacteria', 'Bacteria', 'Terrabacteria group', 'FCB group']:
        return 'Bacteria'
    elif t['group']  in ['Archaea', 'Euryarchaeota', 'TACK group']:
        return 'Archaea'
    elif t['group']  == 'others':
        return others2group[t['tax']]
    else:
        return t['group']

def main():
    tax = parse_placements(glob.glob("placements/*.csv"))
    eur2NOG = parse_hierarchy('COG2eurNOG.tsv')
    taxids = [line.strip() for line in open("Methanotecta.taxid")]
    failed_refs = [line.split()[0] for line in open('references_summary.txt') if 'TreeError' in line]

    tax = include_remaining_queries(tax, eur2NOG, taxids)
    clst2node, node2tab = get_clst2node(glob.glob('tables-gains/*.tab'))
    loss = parse_losses(glob.glob('tables-losses/*.tab'))

    annot = {line.split('\t')[0]:";".join(eval(line.split('\t')[4])) for line in open("eurNOG_annotations.tsv")}
    annot = pd.DataFrame.from_dict(annot, orient="index", columns=["cat"])
    # node = pd.DataFrame.from_dict(clst2node, orient="index", columns=["node"])
    node = pd.DataFrame()
    for k,v in clst2node.items():
        for n in v:
            node = node.append(pd.DataFrame({"node":n}, index=[k]))


    print("manually adding 0M14X with ref 0ZFRM")
    tax = tax.append(pd.DataFrame({"COG":"0ZFRM", "tax":"de-novo genes"}, index=["0M14X"]))


    print_raw_results(tax, annot, node)
    freq = convert_merge_dfs(tax, annot, node, loss)
    freq.to_csv('all_freqs.csv', sep='\t', index=False)

    node_dict = {g[0]:g[1] for g in freq.groupby(['node'])}
    plot_ancestors(node_dict)
    plot_all(node_dict)


if __name__ == '__main__':
    main()

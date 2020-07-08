from collections import defaultdict
from Bio.KEGG import REST
from urllib.error import HTTPError
from Bio import SeqIO
import ete3

header = ['Species', 'protein length', 'arCOG', 'arCOG geneID', 'arCOG description', 'KO', 'KEGG description', 'KEGG pathway', 'eurNOG', 'eggNOG annot', 'eggNOG cat', 'PROKKA annotation', 'PROKKA gene', 'Pfam', 'Pfam description', 'CAZy', 'CAZy description', 'TCDB', 'TCDB description', 'TIGRFAM', 'TIGRFAM description', 'IPR domain', 'IPR description', 'IPR pathway', 'NR hits', 'NR taxonomy', 'NR taxID']

def def_annot():
    return {'Species':'',
           'CAZy':[],
           'CAZy description':set(),
           'NR hits':[],
           'NR taxonomy':"",
           'NR taxID':"",
           'arCOG':[],
           'arCOG geneID':[],
           'arCOG description':[],
           'KO':[],
           'KEGG description':[],
           'KEGG pathway':set(),
           'eurNOG':'',
           'eggNOG annot':'',
           'eggNOG cat':'',
           'TCDB':[],
           'TCDB description':set(),
           'protein length':'',
           'Pfam':[],
           'Pfam description':[],
           'TIGRFAM':[],
           'TIGRFAM description':[],
           'Gene3D':[],
           'Gene3D description':[],
           'SUPERFAMILY':[],
           'SUPERFAMILY description':[],
           'Hamap':[],
           'Hamap description':[],
           'PIRSF':[],
           'PIRSF description':[],
           'SMART':[],
           'SMART description':[],
           'CDD':[],
           'CDD description':[],
           'Coils':[],
           'Coils description':[],
           'SFLD':[],
           'SFLD description':[],
           'IPR domain':set(),
           'IPR description':set(),
           'IPR pathway':set(),
           'PROKKA annotation':'',
           'PROKKA gene':''}

species = ['bin_117', 'bin_91', 'bin_125', 'bin_160', 'bin_172']
species_hika = {'bin_117':'Ca. Hikarchaeum yamanae Bin3',
'bin_125':'Ca. Hikarchaeum sp. Bin4',
'bin_160':'Ca. Hikarchaeum sp. Bin5',
'bin_172':'Ca. Hikarchaeum yamanae Bin2',
'bin_91':'Ca. Hikarchaeum yamanae Bin1'}

annotations = defaultdict(lambda: def_annot())

def parse_eggnog():
    for sp in species:
        for line in open('{}.emapper.annotations'.format(sp)):
            if not line.startswith('#'):
                line = line.strip().split('\t')
                seq = line[0].split('..')[1]
                annotations[seq]['KO'] = line[6].split(',') if line[6] else []
                annotations[seq]['eurNOG'] = line[10].split('|')[0]
                annotations[seq]['eggNOG cat'] = line[11]
                annotations[seq]['eggNOG annot'] = line[12]
                OGs = {o.split('@')[1]:o.split('@')[0] for o in line[9].split(',')}
                if 'arNOG' in OGs:
                    annotations[seq]['arCOG'] = OGs['arNOG']

def parse_arCOG():
    arcog_annotations = {}
    genes = {}
    for line in open('ar14.arCOGdef.tab'):
        line = line.strip().split('\t')
        arcog_annotations[line[0]] = line[3]
        genes[line[0]] = line[2]
    for k, v in annotations.items():
        if v['arCOG'] in arcog_annotations:
            annotations[k]['arCOG description'] = arcog_annotations[v['arCOG']]
            annotations[k]['arCOG geneID'] = genes[v['arCOG']]

def parse_KEGG():
    ko_desc = {}
    ko_path = {}
    for line in open('KEGG_KO_Pathways.csv'):
        line = line.split('\t')
        ko_desc[line[0]] = line[1]
        ko_path[line[0]] = line[2].strip().split(';')

    for k,v in annotations.items():
        if annotations[k]['KO']:
            for ko in annotations[k]['KO']:
                annotations[k]['KEGG description'].append(ko_desc[ko])
                for p in ko_path[ko]:
                    annotations[k]['KEGG pathway'].add(p)

def parse_CAZy():
    cazy_annot = defaultdict(str)
    for line in open("CAZy/fam_activities.txt"):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            if len(line) > 1:
                cazy_annot[line[0]] = line[1].strip()

    cazy_seq = defaultdict(lambda: {sp:set() for sp in species})
    for sp in species:
        for line in open("CAZy/{}.out".format(sp)):
            if not line.startswith('#'):
                line = line.split()
                seq = line[0].split('..')[1]
                cazy = line[2].replace('.hmm', '')
                annotations[seq]['CAZy'].append(cazy)
                annotations[seq]['CAZy description'].add(cazy_annot[cazy.split('_')[0]])

def parse_TCDB():
    tcdb_annot = {}
    for line in open("TCDB/families.tsv"):
        line = line.strip().split('\t')
        tcdb_annot[line[0]] = line[1].replace('</sup>', '').replace('<sup>', '')

    for sp in species:
        for line in open("TCDB/{}.out".format(sp)):
            line = line.split()
            seq = line[0].split('..')[1]
            if len(line[1].split('|')) == 4:
                fam = line[1].split('|')[3]
                superfam = ".".join(fam.split('.')[0:3])
                annotations[seq]['TCDB'].append(fam)
                annotations[seq]['TCDB description'].add(tcdb_annot[superfam])

def parse_interproscan():
    for sp in species:
        for line in open("INTERPROSCAN/{}.faa.tsv".format(sp)):
            line = line.strip().split('\t')
            seq = line[0].split('..')[1]
            analysis = line[3]
            annotations[seq][analysis].append(line[4])
            annotations[seq]["{} description".format(analysis)].append(line[5])
            for i, k in zip([11, 12, 13], ['IPR domain', 'IPR description', 'IPR pathway']):
                if len(line) >= i + 1:
                    annotations[seq][k].add(line[i])

def parse_prokka():
    for sp in species:
        for line in open("../0_genome_bins_sequences/{}.gff".format(sp)):
            line = line.strip().split('\t')
            if len(line) == 9:
                line_dict = {}
                for elem in line[8].split(';'):
                    elem = elem.split('=')
                    line_dict[elem[0]] = elem[1]
                annotations[line_dict['ID']]['PROKKA annotation'] = line_dict['product']
                if 'gene' in line_dict:
                    annotations[line_dict['ID']]['PROKKA gene'] = line_dict['gene']

def ncbi_get_common_ancestor(taxids, ncbi):
    taxids = [t for t in taxids if t]
    if not taxids:
        return '', ''
    t = ncbi.get_topology(taxids, intermediate_nodes=False)
    tax = ncbi.get_taxid_translator([int(t.name)])[int(t.name)]
    return tax, t.name

def parse_diamond():
    ncbi = ete3.ncbi_taxonomy.NCBITaxa()
    seq2taxids = defaultdict(set)
    for sp in species:
        for line in open('DIAMOND_NR/{}.out'.format(sp)): 
            line = line.strip().split('\t')
            sequence = line[0].split('..')[1]
            annotations[sequence]["NR hits"].append("{}:{}".format(line[1], line[3]))
            for t in line[2].split(';'):
                seq2taxids[sequence].add(t)
    for seq, taxids in seq2taxids.items():
        tax, taxid = ncbi_get_common_ancestor(taxids, ncbi)
        annotations[seq]['NR taxonomy'] = tax
        annotations[seq]['NR taxID'] = taxid

def parse_sequences():
    for sp in species:
        for rec in SeqIO.parse("../0_genome_bins_sequences/{}.faa".format(sp), 'fasta'):
            seq = rec.id.split('..')[1]
            annotations[seq]["protein length"] = str(len(rec.seq))
            annotations[seq]["Species"] = species_hika[sp]

def print_annotations():
    with open("all_annotations_Hikarchaeia.csv", 'w') as out:
        print('\t'.join(['Sequence ID'] + header), file=out)
        for k, v in annotations.items():
            cols = [k]
            for h in header:
                if isinstance(v[h], list) or isinstance(v[h], set):
                    cols.append(';'.join(v[h]))
                else:
                    cols.append(v[h])
            print('\t'.join(cols), file=out)

def main():
    parse_eggnog()
    parse_arCOG()
    parse_KEGG()
    parse_CAZy()
    parse_TCDB()
    parse_prokka()
    parse_interproscan()
    parse_sequences()
    parse_diamond()
    print_annotations()

if __name__ == '__main__':
    main()

# def prepare_KEGG_lookup():
#     KOs = set()
#     for sp in species:
#         for line in open('{}.emapper.annotations'.format(sp)):
#             if not line.startswith('#'):
#                 line = line.strip().split('\t')
#                 if line[6]:
#                     for ko in line[6].split(','):
#                         KOs.add(ko)
#
#     general_pathways = ['path:map01110', 'path:map01100', 'path:map01230', 'path:map01120']
#     KO_path = defaultdict(lambda: {'description':"", 'pathways':set()})
#     for ko in KOs:
#         try:
#             KO_path[ko]['description'] = REST.kegg_list(ko).read().strip().split('\t')[1]
#             paths = REST.kegg_link('pathway', ko).read().strip()
#             if paths:
#                 for path in paths.split('\n'):
#                     path = path.split('\t')
#                     if path[1].startswith('path:map') and not any([path[1] == p for p in general_pathways]):
#                         path_str = ':'.join(REST.kegg_list(path[1]).read().strip().split('\t'))
#                         KO_path[ko]['pathways'].add(path_str)
#         except:
#             print("{} did not finish, errors occured".format(ko))
#             KO_path[ko]
#
#     with open('KEGG_KO_Pathways.csv', 'w') as out:
#         for k,v in KO_path.items():
#             print("{}\t{}\t{}".format(k, v['description'], ';'.join(v['pathways'])), file=out)

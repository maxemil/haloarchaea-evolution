# load the eurNOG database so that many genomes can be annotated simoultanously
emapper.py --database eurNOG --servermode

# process all proteomes parallel
nextflow scripts/emapper-parallel.nf --fasta "proteomes/*" --outdir annotations_4.5.1_eurNOG --server_adress eurNOG:localhost:51802

# make clusters based on annotation and  cluster all non-annotated sequences de-novo
python3 scripts/make_EggNOG_clusters.py -f "proteomes/*" -a "annotations_4.5.1_eurNOG/*" -o cluster_eurNOG_4.5.1 -s HiFiX.sif --references

# find and split inconsistent eurNOG (potential gene fusions)
python3 scripts/find_inconsistent_NOGs.py > inconsistent_eurNOGs.txt # with this you can produce data/inconsistent_eurNOGs.txt


mkdir cluster_to_split
for c in $(cut inconsistent_eurNOGs.txt -f1); do mv clusters_faa_4.5.1_eurNOG/$c.fasta clusters_to_split; done

mkdir COG
mkdir HMM_output
python3 scripts/align_COGs.py # to make HMM profiles from the original COG clusters and then run HMMsearch for each cluster with multiple associations to NOGs against the profiles for these


mkdir clusters_splitted
python3 scripts/split_clusters_domain.py -h HMM_output "cluster_to_split/*" # and finally split the clusters based on these searches
## copy all final clusters to cluster_eurNOG_4.5.1_resolved/

# align clusters with 4 or more members
nextflow run scripts/align-clusters.nf --cluster 'cluster_eurNOG_4.5.1_resolved/*.fasta' --output_alignments clusters_eurNOG/alignments --output_faa clusters_eurNOG/faa -qs 20

# generate trees and ufboot samples for all alignemnts
nextflow run scripts/batch-run-trees.nf --alignments 'clusters_eurNOG/alignments/*' --output_trees clusters_eurNOG/trees -resume

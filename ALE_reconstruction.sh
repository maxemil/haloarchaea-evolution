# clone the nextflow pipeline for running ALE
git clone https://github.com/maxemil/ALE-pipeline.git

# build singularity image of ALE
sudo singularity build ALE.sif ALE-pipeline/Singularity

# create species name abbreviations to use in ALE
python3 ALE-pipeline/bin/create_name_codes.py -i species_names.txt -o species_map.txt

# make sure the names in the species tree are the same as in the clusters, 
# otherwise create another mapping file for those, genes_map.txt

# the file fraction_missing.txt should contain completeness estimates,
# e.g 0.2 for 80% complete

# run ALE using ufboots from clusters and small clusters 
nextflow ALE-pipeline/main.nf \
        --species_tree_files species_trees_input/species.tree \
        --input_files "clusters_eurNOG/trees/" \
        --input_extension ".ufboot" \
        --small_cluster "clusters_eurNOG/faa/cluster_small/*" \
        --genes_map genes_map.txt \
        --species_map species_map.txt \
        --fraction_missing fraction_missing.txt \
        --output_ale MGIV-eurNOG \
        --output_samples ALE-samples \
        --separator ".." \
        --outgroup_taxa '["Euryarchaeota_archaeon_JdFR_21","Archaeoglobus_fulgidus_DSM_4304","Archaeoglobi_archaeon_JdFR_42"]' \
        --python3 /usr/bin/python \
        -resume -qs 8

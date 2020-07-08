git clone https://github.com/maxemil/originations-placement.git

sudo singularity build originations-placement.sif originations-placement/Singularity

mkdir tables-gains
cp MGIV-eurNOG-tables/tab_files/gains* tables-gains


# download NOG2EurNOG.tab

mkdir references_trimal
mkdir failed_references_trimal

nextflow run originations-placement/main.nf \
        --fasta 'queries_split_all/*' \
        --reference_alignments 'references_trimal/*.trimmed_alg' \
        --reference_trees 'references_trimal/*.tree' \
        --precompiled_taxa precompiled.tax.map \
        --og_hierarchy COG2eurNOG_split.tsv \
        --remove_taxids Methanotecta.taxid \
        --placment_dir placements_addfragments \
        --references_trimmed_dir references_trimmed_trimal \
        --queries_dir halo_queries_addfragments \
        -with-singularity originations-placement.sif \
        -resume

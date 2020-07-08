mkdir CAZy
mkdir TCDB
mkdir INTERPROSCAN
mkdir DIAMOND_NR

function annotate {
  f=$1
  base=$(basename ${f%%.faa})
  hmmsearch --tblout CAZy/$base.out -E 1e-3 CAZy/dbCAN-HMMdb-V8.txt $f
  blastp -num_threads 4 -outfmt 6 \
            -query $f -db TCDB/TCDB.fasta \
            -out TCDB/$base.out -evalue 1e-20
  blastp -num_threads 4 -outfmt 6 \
            -query $f -db VFDB/VFDB_setB_pro.fas \
            -out VFDB/$base.out -evalue 1e-20
  /local/two/Software/interproscan/interproscan-5.42-78.0/interproscan.sh -i $f \
            -f tsv \
            -t p \
            -appl TIGRFAM,SMART,SUPERFAMILY,Pfam,Gene3D,Hamap,SFLD,CDD,Coils,PIRSF \
            --pathways \
            --iprlookup \
            --output-dir INTERPROSCAN \
            --cpu 4
  # diamond blastp -q $faa --db $DIAMOND_NR -e 1e-5 --out DIAMOND/"$faa".out \
  #           --outfmt 6 qseqid sseqid staxids evalue pident length mismatch gapopen qstart qend sstart send bitscore \
  #           --threads 15 --top 2
}


ls proteomes/*.faa  | parallel -j2 annotate

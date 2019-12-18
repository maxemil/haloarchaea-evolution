docker run --rm -v $PWD:$PWD \
      maxemil/phyloskeleton phyloSkeleton.pl \
      -o . \
      -p eury \
      -H data/archaea56.hmm \
      --best-match-only \
      -c 0.6 \
      -C 20 \
      -s data/eury.selection.user_refined.tab \
      -t eury.taxids.tab
      # -s eury_bins.tab \
      # -n Archaea.csv \
      # -J Euryarchaeota_JGI.tab \
      # -l selLevels.tab \

cp ../1_taxon_sampling_Eury/eury.fasta/arCOG*.fasta arCOG_extracted

for aln in arCOG_extracted/*.fasta;
do
  mafft-einsi --thread 10 $aln > ${aln%.*}".aln";
  trimal -in  ${aln%.*}".aln" -out ${aln%.*}".trim" -fasta -gappyout;
done

concatenateRenameAlignment.pl \
                  -m data/eury.map \
                  -s .trim \
                  arCOG_extracted/*.trim > euryarchaeota_TACK_56_arCOGs.fasta

iqtree \
  -s euryarchaeota_TACK_56_arCOGs.fasta \
  -m LG+I+G4 \
  -nt AUTO \
  -pre euryarchaeota_TACK \
  -fast
  # -bb 1000
  # -m TEST \
  # -redo \
  # -mredo


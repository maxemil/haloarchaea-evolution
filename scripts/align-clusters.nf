params.cluster = ""
params.output_alignments = "alignments"
params.output_faa = "cluster_faa"
params.output_trees = "trees"

cluster = Channel.fromPath(params.cluster)


params.alignments = ""
mafft_alignments = Channel.fromPath(params.alignments)

// process sortClustersBySize {
//   input:
//   file faa from cluster
//
//   output:
//   file "small_cluster/${faa.simpleName}.faa" into small_faa optional true
//   file "${faa.simpleName}.faa" into normal_faa optional true
//
//   publishDir "${params.output_faa}", mode: 'copy', pattern: 'small_cluster/*'
//   tag {"${faa.simpleName}"}
//
//   script:
//   """
//   if [ "\$(grep -c '^>' $faa)" -lt 4 ]
//   then
//     mkdir small_cluster
//     mv $faa small_cluster/${faa.simpleName}.faa
//   else
//     cp -L $faa ${faa.simpleName}.faa
//   fi
//   """
// }

cluster.into{normal_faa_alignment; normal_faa_filering}

// process prequal {
//   input:
//   file faa from normal_faa_alignment
//
//   output:
//   file "${faa.simpleName}.clean.filtered" into prequal_results
//
//   tag {"${faa.simpleName}"}
//   cpus 1
//   publishDir "alignments", mode: 'copy'
//
//   script:
//   """
//   sed 's/*//g' $faa > ${faa.simpleName}.clean
//   prequal ${faa.simpleName}.clean
//   """
// }
//
// process alignSequences {
//   input:
//   file faa from prequal_results
//
//   output:
//   file "${faa.simpleName}.aln" into mafft_alignments
//
//   tag {"${faa.simpleName}"}
//   cpus 6
//   publishDir "alignments", mode: 'copy'
//
//   script:
//   """
//   mafft-einsi --anysymbol --thread ${task.cpus} $faa > ${faa.simpleName}.aln
//   """
// }


filter_alignments = normal_faa_filering.combine(mafft_alignments)

process filterAlignmentSize{
  input:
  set file(faa), file(aln) from filter_alignments

  output:
  file "$aln" into filtered_divvyied_alignments optional true
  file "$faa" into filtered_faa optional true
  file "small_cluster/*.faa" into filtered_small_cluster optional true

  publishDir "${params.output_alignments}", mode: 'copy', pattern: '*.aln'
  publishDir "${params.output_faa}", mode: 'copy', pattern: 'small_cluster/*'
  publishDir "${params.output_faa}", mode: 'copy', pattern: '*.fasta'
  stageInMode 'copy'
  tag {"${faa.simpleName}"}

  when:
  "${faa.simpleName}" == "${aln.simpleName}"

  script:
  """
  #! /usr/bin/env python3
  from Bio import SeqIO
  import os

  os.makedirs("small_cluster")
  missing_data = ['-', 'X']

  aln = {rec.id:rec for rec in SeqIO.parse("$aln", 'fasta')}
  faa = {rec.id:rec for rec in SeqIO.parse("$faa", 'fasta')}
  all_missing_data = []
  for recid, rec in aln.items():
        if all([c in missing_data for c in rec.seq]):
            all_missing_data.append(recid)
  [aln.pop(recid) for recid in all_missing_data]
  small = "" if len(aln.keys()) > 3 else 'small_cluster/'
  with open("{}$faa".format(small), 'w') as faa_file, open("{}$aln".format(small), 'w') as divvy_file:
      for recid, rec in faa.items():
          if recid in aln.keys():
              SeqIO.write(aln[recid], divvy_file, 'fasta')
              SeqIO.write(rec, faa_file, 'fasta')
          else:
              with open("small_cluster/{}.faa".format(recid.replace('.', '_')), 'w') as single:
                  SeqIO.write(rec, single, 'fasta')
  if small:
    os.remove("$aln")
    os.remove("$faa")
  """

}


process cluster2Tree {
  input:
  file aln from filtered_divvyied_alignments

  output:
  file "${aln.simpleName}.*" into treefiles

  publishDir "${params.output_trees}", mode: 'copy'
  tag {"${aln.simpleName}"}
  cpus 4

  script:
  """
  iqtree -s $aln \
	 -pre ${aln.simpleName} \
	 -bb 1000 -bnni \
	 -nt ${task.cpus} \
	 -m TESTNEW \
	 -mset LG \
	 -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
	 -keep-ident \
	 -seed 12345 \
	 -wbtl
  """
}

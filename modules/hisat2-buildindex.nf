process buildindexHISAT2{
  label 'mg00_hisat2_buildindex'
  conda params.buildindexHISAT2.conda
  cpus params.resources.buildindexHISAT2.cpus
  memory params.resources.buildindexHISAT2.mem
  queue params.resources.buildindexHISAT2.queue
  clusterOptions params.resources.buildindexHISAT2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_hisat2_buildindex", mode: 'symlink'
  input:
  val(index_dir)
  val(index_name)
  path(fasta)
  path(annot)
  
  
  output:
  tuple(env(index_base), path(index_dir))

  shell:
  '''
  currdir=$(pwd)
  index_base=$currdir/!{index_dir}/!{index_name}
  mkdir !{index_dir}

  hisat2_extract_splice_sites.py !{annot} > genome.ss
  hisat2_extract_exons.py !{annot} > genome.exon 

  hisat2-build -p !{params.resources.buildindexHISAT2.cpus} \
        !{fasta} \
        !{params.buildindexHISAT2.options} \
        --ss genome.ss --exon genome.exon $index_base
        
  '''
}
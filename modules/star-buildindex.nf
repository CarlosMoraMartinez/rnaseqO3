process buildindexSTAR{
  label 'mg00_star_buildindex'
  conda params.buildindexSTAR.conda
  cpus params.resources.buildindexSTAR.cpus
  memory params.resources.buildindexSTAR.mem
  queue params.resources.buildindexSTAR.queue
  clusterOptions params.resources.buildindexSTAR.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_star_buildindex", mode: 'symlink'
  input:
  val(index_dir)
  path(fasta)
  path(annot)
  
  output:
  path(index_dir)

  shell:
  '''
  mkdir !{index_dir}
  STAR --runThreadN !{params.resources.buildindexSTAR.cpus} \
--runMode genomeGenerate \
--genomeDir !{index_dir} \
--genomeFastaFiles !{fasta} \
--sjdbGTFfile !{annot} \
--genomeSAindexNbases !{params.buildindexSTAR.genomeSAindexNbases} \
--sjdbOverhang !{params.buildindexSTAR.sjdbOverhang} #99 o 100
  '''
}
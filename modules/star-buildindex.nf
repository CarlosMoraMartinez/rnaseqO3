process buildindexSTAR{
  label 'mg02_star_buildindex'
  conda params.buildindexSTAR.conda
  cpus params.resources.buildindexSTAR.cpus
  memory params.resources.buildindexSTAR.mem
  queue params.resources.buildindexSTAR.queue
  clusterOptions params.resources.buildindexSTAR.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_star_buildindex", mode: 'symlink'
  input:
  tuple(index_name)
  
  output:
  path(index_name)

  shell:
  '''
  mkdir !{index_name}
  STAR --runThreadN !{params.resources.buildindexSTAR.cpus} \
--runMode genomeGenerate \
--genomeDir !{index_name} \
--genomeFastaFiles !{params.resources.alignSTAR.genome_fasta} \
--sjdbGTFfile !{params.resources.alignSTAR.genome_gtf} \
--genomeSAindexNbases !{params.resources.alignSTAR.genomeSAindexNbases} \
--sjdbOverhang !{params.resources.alignSTAR.sjdbOverhang} #99 o 100
  '''
}
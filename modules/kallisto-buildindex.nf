process buildindexKallisto{
  label 'mg00_kallisto_buildindex'
  conda params.buildindexKallisto.conda
  cpus params.resources.buildindexKallisto.cpus
  memory params.resources.buildindexKallisto.mem
  queue params.resources.buildindexKallisto.queue
  clusterOptions params.resources.buildindexKallisto.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_kallisto_buildindex", mode: 'symlink'
  input:
  val(index_name)
  path(transcripts_fasta)
  
  output:
  path("*.idx")

  shell:
  '''
  kallisto index -i !{index_name}.idx \
    --threads !{params.resources.buildindexKallisto.cpus} \
    -k !{params.buildindexKallisto.k} !{params.buildindexKallisto.options} !{transcripts_fasta}
  '''
}
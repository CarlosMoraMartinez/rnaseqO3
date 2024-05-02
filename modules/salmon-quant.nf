process quantSalmon{
  label 'mg03_salmon_quant'
  conda params.quantSalmon.conda
  cpus params.resources.quantSalmon.cpus
  memory params.resources.quantSalmon.mem
  queue params.resources.quantSalmon.queue
  clusterOptions params.resources.quantSalmon.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_salmon_quant", mode: 'symlink'
  input:
  path(salmon_index)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*_salmon'))

  shell:
  '''
  outdir=!{sample_id}_salmon
  errorfile=!{sample_id}.hisat2.err

  salmon quant --threads !{params.resources.quantSalmon.cpus} \
      -i !{salmon_index} -l IU !{params.quantSalmon.options} \
      -1 !{fastq[0]} -2 !{fastq[1]} \
      --validateMappings -o $outdir 2> $errorfile

  '''
}
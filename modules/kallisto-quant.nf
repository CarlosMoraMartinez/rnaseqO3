process quantKallisto{
  label 'mg03_kallisto_quant'
  conda params.quantKallisto.conda
  cpus params.resources.quantKallisto.cpus
  memory params.resources.quantKallisto.mem
  queue params.resources.quantKallisto.queue
  clusterOptions params.resources.quantKallisto.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_kallisto_quant", mode: 'symlink'
  input:
  path(kallisto_index)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*_kallisto'), path("*.kallisto.err"))

  shell:
  '''
  outdir=!{sample_id}_kallisto
  errorfile=!{sample_id}.kallisto.err

  kallisto quant -i !{kallisto_index} \
    --threads !{params.resources.quantKallisto.cpus} \
    -b !{params.quantKallisto.num_bootstraps} !{params.quantKallisto.options} \
    -o $outdir \
    !{fastq[0]} !{fastq[1]} 2> $errorfile

  '''
}
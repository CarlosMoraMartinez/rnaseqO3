process quantBamSalmon{
  label 'mg03_salmon_quant_frombam'
  conda params.quantBamSalmon.conda
  cpus params.resources.quantBamSalmon.cpus
  memory params.resources.quantBamSalmon.mem
  queue params.resources.quantBamSalmon.queue
  clusterOptions params.resources.quantBamSalmon.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_salmon_quant_frombam", mode: 'symlink'
  input:
  tuple(val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val("saln-${flag}"), val(sample_id), path('*_salmonaln'))

  shell:
  '''
  outdir=!{sample_id}_!{flag}_salmonaln
  errorfile=!{sample_id}_!{flag}.salmonaln.err

  ## Quantifying in alignment based mode
  salmon quant --threads !{params.resources.quantBamSalmon.cpus} \
      -t !{params.quantBamSalmon.transcripts_fasta} \
      -l IU !{params.quantBamSalmon.options} \
      -a !{bam} \
      -o $outdir 2> $errorfile

  '''
}
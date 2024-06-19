process mergeSalmon{
  label 'mg03_salmon_merge'
  conda params.mergeSalmon.conda
  cpus params.resources.mergeSalmon.cpus
  memory params.resources.mergeSalmon.mem
  queue params.resources.mergeSalmon.queue
  clusterOptions params.resources.mergeSalmon.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_salmon_merge", mode: 'symlink'
  input:
  tuple(val(tag), val(sample_names), path(salmon_dirs))
 
  output:
  tuple(path('*_numreads.tsv'), path('*_tpm.tsv'))

  shell:
  '''
  salmon quantmerge --quants !{salmon_dirs} --names !{sample_names} --column NumReads --output !{tag}_numreads.tsv
  salmon quantmerge --quants !{salmon_dirs} --names  !{sample_names} --column TPM --output !{tag}_tpm.tsv
  '''
}
process buildindexSubread{
  label 'mg00_subread_buildindex'
  conda params.buildindexSubread.conda
  cpus params.resources.buildindexSubread.cpus
  memory params.resources.buildindexSubread.mem
  queue params.resources.buildindexSubread.queue
  clusterOptions params.resources.buildindexSubread.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_subread_buildindex", mode: 'symlink'
  input:
  val(index_dir)
  val(index_name)
  path(fasta)
  
  output:
  tuple(env(index_base), path(index_dir))

  shell:
  '''
  currdir=$(pwd)
  index_base=$currdir/!{index_dir}/!{index_name}
  mkdir !{index_dir}

  subread-buildindex !{fasta} !{params.buildindexSubread.options} -o !{index_dir}/!{index_name} 2> !{index_dir}/!{index_name}'.err'
  '''
}
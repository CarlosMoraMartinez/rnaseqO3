process buildindexBBMap{
  label 'mg00_bbmap_buildindex'
  conda params.buildindexBBMap.conda
  cpus params.resources.buildindexBBMap.cpus
  memory params.resources.buildindexBBMap.mem
  queue params.resources.buildindexBBMap.queue
  clusterOptions params.resources.buildindexBBMap.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_bbmap_buildindex", mode: 'symlink'
  input:
  path(fasta)
  
  output:
  tuple(path("ref"))

  shell:
  '''
  mem=$(echo !{params.resources.buildindexBBMap.mem} | sed "s/ GB/g/")

  sed !{params.buildindexBBMap.clean_fasta_regex} !{params.buildindexBBMap.fasta} >tmp.fasta
  bbmap.sh -Xmx$mem \
    ref=tmp.fasta \
    k=!{params.buildindexBBMap.k} 2> bbmap_index.err

  grep ">" tmp.fasta > new_seq_names.txt # to check that reference names are ok
  rm tmp.fasta

  '''
}
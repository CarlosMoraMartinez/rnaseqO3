
process buildRefflatFromGTF{
  label 'mg00_build_refflat'
  conda params.buildRefflatFromGTF.conda
  cpus params.resources.buildRefflatFromGTF.cpus
  memory params.resources.buildRefflatFromGTF.mem
  queue params.resources.buildRefflatFromGTF.queue
  clusterOptions params.resources.buildRefflatFromGTF.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_build_refflat", mode: 'symlink'
  input:
  path(genome_fasta)
  path(annot)
  
  output:
  path('*.refflat')

  shell:
  '''
  outname=$(basename -s .gtf !{annot})'.refflat'
        
  java -jar !{params.software.gtftorefflat_path} \
  -g !{annot} \
  -R !{genome_fasta} \
  -r $outname 2> gtftorefflat.err

  '''
}


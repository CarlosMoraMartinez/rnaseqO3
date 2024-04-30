process mapdecoySalmontools{
  label 'mg00_salmontools_mapdecoy'
  conda params.mapdecoySalmontools.conda
  cpus params.resources.mapdecoySalmontools.cpus
  memory params.resources.mapdecoySalmontools.mem
  queue params.resources.mapdecoySalmontools.queue
  clusterOptions params.resources.mapdecoySalmontools.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_salmontools_mapdecoy", mode: 'symlink'
  input:
  val(name_decoy)
  path(genome_fasta)
  path(transcripts_fasta)
  path(annot)
  
  output:
  tuple(path("${name_decoy}/gentrome.fa"), path("${name_decoy}/decoys.txt"))

  shell:
  '''
  generateDecoyTranscriptome.sh -j !{params.resources.mapdecoySalmontools.cpus} -a !{annot} -g !{genome_fasta} -t !{transcripts_fasta} -o name_decoy
 
  '''
}
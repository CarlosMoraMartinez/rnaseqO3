
process buildIntervalListFromBed{
  label 'mg00_build_intervallist'
  conda params.buildIntervalListFromBed.conda
  cpus params.resources.buildIntervalListFromBed.cpus
  memory params.resources.buildIntervalListFromBed.mem
  queue params.resources.buildIntervalListFromBed.queue
  clusterOptions params.resources.buildIntervalListFromBed.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_build_intervallist", mode: 'symlink'
  input:
  path(genome_fasta)
  path(bed)
  
  output:
  path('*.interval_list')

  shell:
  '''
  dict_file=$(basename -s '.fa' !{genome_fasta})'.dict'
  outfile=$(basename -s '.bed' !{bed})'.interval_list'

  # First create a .dict file from the genome
  java -jar !{params.software.gatk_path} CreateSequenceDictionary -R !{genome_fasta} 2> CreateSequenceDictionary.err

  java -jar !{params.software.picard_path} BedToIntervalList \
    I=!{bed} \
    O=$outfile \
    SD=$dict_file 2> BedToIntervalList.err
        
  '''
}


process doTrimmomatic{
  label 'mg01_trimmomatic'
  conda params.doTrimmomatic.conda
  cpus params.resources.doTrimmomatic.cpus
  memory params.resources.doTrimmomatic.mem
  queue params.resources.doTrimmomatic.queue
  clusterOptions params.resources.doTrimmomatic.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg01_trimmomatic", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(fastq))
  
  output:
  tuple(val(illumina_id), path('*.trimmtc.fastq.gz'), path('*.trimmtc.single.fastq.gz'), path('*.trimlog'), path('*.trimlog.out'))

  shell:
  '''
  trimmomatic PE -threads !{params.resources.doTrimmomatic.cpus} -phred33 \
        -trimlog !{illumina_id}.trimlog \
        !{fastq[0]} !{fastq[1]} \
        !{illumina_id}_1.trimmtc.fastq.gz \
        !{illumina_id}_1.trimmtc.single.fastq.gz \
        !{illumina_id}_2.trimmtc.fastq.gz \
        !{illumina_id}_2.trimmtc.single.fastq.gz \
        ILLUMINACLIP:!{params.doTrimmomatic.illuminaclip} \
        SLIDINGWINDOW:!{params.doTrimmomatic.slidingwindow} \
        MINLEN:!{params.doTrimmomatic.minlen} 2> !{illumina_id}.trimlog.out
  '''
}
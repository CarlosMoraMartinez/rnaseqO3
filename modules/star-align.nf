process alignSTAR{
  label 'mg02_star_align'
  conda params.alignSTAR.conda
  cpus params.resources.alignSTAR.cpus
  memory params.resources.alignSTAR.mem
  queue params.resources.alignSTAR.queue
  clusterOptions params.resources.alignSTAR.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_star_align", mode: 'symlink'
  input:
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*.bam'), path('*Log.final.out'), path('*Log.out'), path('*Log.progress.out'))

  shell:
  '''

  STAR --runThreadN !{params.resources.alignSTAR.cpus} !{params.alignSTAR.options} \
  --readFilesCommand zcat \
  --genomeDir !{params.resources.alignSTAR.index} \
  --readFilesIn !{fastq[0]},!{fastq[1]}  \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix !{sample_id}_
  samtools index !{sample_id}_Aligned.sortedByCoord.out.bam

  '''
  stub:
  """
  touch $sample_id'_Aligned.sortedByCoord.out.bam'
  touch $sample_id'_Log.final.out'
  touch $sample_id'_Log.out'
  touch $sample_id'_Log.progress.out'
  """
}
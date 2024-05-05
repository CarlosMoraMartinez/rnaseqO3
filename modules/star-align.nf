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
  path(star_index_dir)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), 
        path('*Aligned.sortedByCoord.out.bam'), 
        path('*Aligned.sortedByCoord.out.bam.bai'), 
        path('*Aligned.toTranscriptome.sorted.bam'), 
        path('*Aligned.toTranscriptome.sorted.bam.bai'), 
        path("*_SJ.out.tab"), 
        path('*Log.final.out'), 
        path('*Log.out'), 
        path('*Log.progress.out'),
        path('*ReadsPerGene.out.tab'),
        path('*star.flagstat'),
        path('*star.transcriptome.flagstat'))

  shell:
  '''

  STAR --runThreadN !{params.resources.alignSTAR.cpus} !{params.alignSTAR.options} \
  --quantMode TranscriptomeSAM GeneCounts \
  --readFilesCommand zcat \
  --genomeDir !{star_index_dir} \
  --readFilesIn !{fastq[0]} !{fastq[1]}  \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix !{sample_id}_
  
  samtools index !{sample_id}_Aligned.sortedByCoord.out.bam

  sortedbamtranscr=!{sample_id}_Aligned.toTranscriptome.sorted.bam
  samtools sort !{sample_id}_Aligned.toTranscriptome.out.bam -O bam -o $sortedbamtranscr
  samtools index $sortedbamtranscr
  rm !{sample_id}_Aligned.toTranscriptome.out.bam

  flagstat=!{sample_id}.star.flagstat
  samtools flagstat !{sample_id}_Aligned.sortedByCoord.out.bam > $flagstat

  flagstat_trans=!{sample_id}.star.transcriptome.flagstat
  samtools flagstat $sortedbamtranscr > $flagstat_trans

  '''
  stub:
  """
  touch $sample_id'_Aligned.sortedByCoord.out.bam'
  touch $sample_id'__SJ.out.tab'
  touch $sample_id'_Log.final.out'
  touch $sample_id'_Log.out'
  touch $sample_id'_Log.progress.out'
  touch $sample_id'.star.flagstat'
  """
}
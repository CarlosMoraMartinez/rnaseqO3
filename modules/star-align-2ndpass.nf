process alignSTAR2ndPass{
  label 'mg02_star_align_2nd_pass'
  conda params.alignSTAR2ndPass.conda
  cpus params.resources.alignSTAR2ndPass.cpus
  memory params.resources.alignSTAR2ndPass.mem
  queue params.resources.alignSTAR2ndPass.queue
  clusterOptions params.resources.alignSTAR2ndPass.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_star_align_2nd_pass", mode: 'symlink'
  input:
  path(star_index_dir)
  tuple(val(sample_id), path(fastq))
  path(splice_junctions_all)
  
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
        path('*starp2.flagstat'),
        path('*starp2.transcriptome.flagstat'))

  shell:
  '''
  
  STAR --runThreadN !{params.resources.alignSTAR2ndPass.cpus} !{params.alignSTAR2ndPass.options} \
  --readFilesCommand zcat \
  --quantMode TranscriptomeSAM GeneCounts \
  --genomeDir !{star_index_dir} \
  --readFilesIn !{fastq[0]},!{fastq[1]}  \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix !{sample_id}_p2_ \
  --sjdbFileChrStartEnd !{splice_junctions_all}

  samtools index !{sample_id}_p2_Aligned.sortedByCoord.out.bam

  sortedbamtranscr=!{sample_id}_p2_Aligned.toTranscriptome.sorted.bam
  samtools sort !{sample_id}_p2_Aligned.toTranscriptome.out.bam -O bam -o $sortedbamtranscr
  samtools index $sortedbamtranscr
  rm !{sample_id}_p2_Aligned.toTranscriptome.out.bam

  flagstat=!{sample_id}.starp2.flagstat
  samtools flagstat !{sample_id}_p2_Aligned.sortedByCoord.out.bam > $flagstat

  flagstat_trans=!{sample_id}.starp2.transcriptome.flagstat
  samtools flagstat $sortedbamtranscr > $flagstat_trans

  '''
  stub:
  """
  touch $sample_id'_Aligned.sortedByCoord.out.bam'
  touch $sample_id'_Aligned.sortedByCoord.out.bam.bai'
  touch $sample_id'_Aligned.toTranscriptome.sorted.bam'
  touch $sample_id'_Aligned.toTranscriptome.sorted.bam.bai'
  touch $sample_id'__SJ.out.tab'
  touch $sample_id'_Log.final.out'
  touch $sample_id'_Log.out'
  touch $sample_id'_Log.progress.out'
  touch $sample_id'.starp2.flagstat'
  touch $sample_id'.starp2.transcriptome.flagstat'
  """
}
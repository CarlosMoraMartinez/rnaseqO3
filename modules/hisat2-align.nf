process alignHISAT2{
  label 'mg02_hisat2_align'
  conda params.alignHISAT2.conda
  cpus params.resources.alignHISAT2.cpus
  memory params.resources.alignHISAT2.mem
  queue params.resources.alignHISAT2.queue
  clusterOptions params.resources.alignHISAT2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_hisat2_align", mode: 'symlink'
  input:
  val(hisat2_index_base)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*.sorted.bam'), path('*.bai'), path('*.hisat2.flagstat'), path('*.hisat2.summary.txt'), path('*.hisat2.metfile.txt'), path('*.hisat2.newsplices.txt'))

  shell:
  '''
  outbam=!{sample_id}.hisat2.sorted.bam
  flagstat=!{sample_id}.hisat2.flagstat
  summaryfile=!{sample_id}.hisat2.summary.txt
  metrics_file=!{sample_id}.hisat2.metfile.txt
  newsplicesites_file=!{sample_id}.hisat2.newsplices.txt
  errorfile=!{sample_id}.hisat2.err

  hisat2 --threads !{params.resources.alignHISAT2.cpus} \
        !{params.alignHISAT2.options} \
        --new-summary --summary-file $summaryfile \
        --novel-splicesite-outfile $newsplicesites_file \
        --met-file $metrics_file \
        -x !{hisat2_index_base}  \
        -1 !{fastq[0]} -2 !{fastq[1]}  2>$errorfile |  \
        tee >(samtools flagstat - > $flagstat) |  \
        samtools sort -O BAM -o $outbam
  samtools index $outbam

  '''
}
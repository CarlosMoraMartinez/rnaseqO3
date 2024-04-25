process alignHISAT2{
  label 'mg02_hisat2align'
  conda params.alignHISAT2.conda
  cpus params.resources.alignHISAT2.cpus
  memory params.resources.alignHISAT2.mem
  queue params.resources.alignHISAT2.queue
  clusterOptions params.resources.alignHISAT2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_hisat2align", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(fastq))
  
  output:
  tuple(val(illumina_id), path('*.sorted.bam'), path('*.hisat2.flagstat'), path('*.hisat2.summary.txt'), path('*.hisat2.metfile.txt'), path('*.hisat2.newsplices.txt'))

  shell:
  '''
  outbam=!{illumina_id}.hisat2.sorted.bam
  flagstat=!{illumina_id}.hisat2.flagstat
  summaryfile=!{illumina_id}.hisat2.summary.txt
  metrics_file=!{illumina_id}.hisat2.metfile.txt
  newsplicesites_file=!{illumina_id}.hisat2.newsplices.txt
  errorfile=!{illumina_id}.hisat2.err

  hisat2 --threads !{params.resources.alignHISAT2.cpus} \
        !{params.alignHISAT2.options} \
        --new-summary --summary-file $summaryfile \
        --novel-splicesite-outfile $newsplicesites_file \
        --met-file $metrics_file \
        -x /home/ccarlos/data/resources/indices/hisat2_mm38_Kim_genome_trans/grcm38_tran/genome_tran  \
        -1 !{fastq[0]} -2 !{fastq[1]}  2>$errorfile |  \
        tee >(samtools flagstat - > $flagstat) |  \
        samtools sort -O BAM -o $outbam

  '''
}
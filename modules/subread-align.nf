process alignSubread{
  label 'ubread_align'
  conda params.alignSubread.conda
  cpus params.resources.alignSubread.cpus
  memory params.resources.alignSubread.mem
  queue params.resources.alignSubread.queue
  clusterOptions params.resources.alignSubread.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_subread_align", mode: 'symlink'
  input:
  val(subread_index_prefix)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*.sorted.bam'), path('*.bai'), path('*.indel.vcf'), path('*.subread.flagstat'), path('*.subread.err'))

  shell:
  '''
  outbam=!{sample_id}.subread.sorted.bam
  flagstat=!{sample_id}.subread.flagstat
  errorfile=!{sample_id}.subread.err

  subread-align -T !{params.resources.alignSubread.cpus} \
  -t 0 !{params.alignSubread.options} \
  -a !{params.alignSubread.annot} --sortReadsByCoordinates \
  -i !{subread_index_prefix} \
  -r !{fastq[0]} -R !{fastq[1]} -o $outbam 2>$errorfile

  samtools flagstat $outbam > $flagstat
  '''
}
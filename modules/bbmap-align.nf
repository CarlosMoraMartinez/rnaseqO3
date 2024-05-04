process alignBBMap{
  label 'ubread_align'
  conda params.alignBBMap.conda
  cpus params.resources.alignBBMap.cpus
  memory params.resources.alignBBMap.mem
  queue params.resources.alignBBMap.queue
  clusterOptions params.resources.alignBBMap.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_bbmap_align", mode: 'symlink'
  input:
  path(bbmap_index_dir)
  tuple(val(sample_id), path(fastq))
  
  output:
  tuple(val(sample_id), path('*.sorted.bam'), path('*.sorted.bam.bai'), path('*.bbmap.flagstat'), path('*.bbmap.err'))

  shell:
  '''
  outbam=!{sample_id}.bbmap.bam
  outbam_sorted=!{sample_id}.bbmap.sorted.bam 
  flagstat=!{sample_id}.bbmap.flagstat
  errorfile=!{sample_id}.bbmap.err

  bbmap.sh in=!{fastq[0]} in2=!{fastq[1]} \
    unpigz=t \
    out=$outbam !{params.alignBBMap.options} \
    maxindel=!{params.alignBBMap.maxindel} \
    secondary=!{params.alignBBMap.secondary}
    ambig=!{params.alignBBMap.ambig} \
    intronlen=!{params.alignBBMap.intronlen}\
    xstag=!{params.alignBBMap.xstag} 2>$errorfile

  samtools sort $outbam -o $outbam_sorted
  samtools index $outbam_sorted

  samtools flagstat $outbam_sorted > $flagstat

  rm $outbam
  '''
}
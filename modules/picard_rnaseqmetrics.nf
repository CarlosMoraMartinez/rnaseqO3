process picardRNASeqMetrics{
  label 'mg04_picard_rnaseqmetrics'
  conda params.picardRNASeqMetrics.conda
  cpus params.resources.picardRNASeqMetrics.cpus
  memory params.resources.picardRNASeqMetrics.mem
  queue params.resources.picardRNASeqMetrics.queue
  clusterOptions params.resources.picardRNASeqMetrics.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg04_picard_rnaseqmetrics", mode: 'symlink'
  input:
  path(refflat)
  path(rRNA_intervallist)
  tuple(val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val(flag), val(sample_id), path('*_picard_rnaseqmetrics.txt'), path('*picard_rnaseqmetrics.err'))

  shell:
  '''
  outreport=!{sample_id}-!{flag}_picard_rnaseqmetrics.txt
  errorfile=!{sample_id}-!{flag}_picard_rnaseqmetrics.err
  newbamname=!{sample_id}_!{flag}.bam
  newbainame=!{sample_id}_!{flag}.bam.bai

  mv !{bam} $newbamname
  mv !{bai} $newbainame
  java -jar !{params.software.picard_path} CollectRnaSeqMetrics \
              I=$newbamname O=$outreport \
              REF_FLAT=!{refflat} !{params.picardRNASeqMetrics.options} \
              RIBOSOMAL_INTERVALS=!{rRNA_intervallist} \
              STRAND=!{params.picardRNASeqMetrics.strand} 2> $errorfile

  '''
}
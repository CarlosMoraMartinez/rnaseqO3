process rnaSeQC{
  label 'mg04_rnaseqc'
  conda params.rnaSeQC.conda
  cpus params.resources.rnaSeQC.cpus
  memory params.resources.rnaSeQC.mem
  queue params.resources.rnaSeQC.queue
  clusterOptions params.resources.rnaSeQC.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg04_rnaseqc", mode: 'symlink'
  input:
  path(annot)
  tuple(val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val(flag), val(sample_id), path('*_rnaseqc'), path('*_rnaseqc.err'))

  shell:
  '''
  outdir=!{sample_id}_!{flag}_rnaseqc
  errorfile=!{sample_id}_!{flag}_rnaseqc.err

  mkdir $outdir

  rnaseqc !{params.rnaSeQC.options} \
   --stranded !{params.rnaSeQC.stranded} \
   --detection-threshold !{params.rnaSeQC.detection_threshold} \
   !{annot} !{bam} $outdir 2> $errorfile 

  '''
}
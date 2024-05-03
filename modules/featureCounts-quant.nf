process quantFeatureCounts{
  label 'mg03_FeatureCounts_quant'
  conda params.quantFeatureCounts.conda
  cpus params.resources.quantFeatureCounts.cpus
  memory params.resources.quantFeatureCounts.mem
  queue params.resources.quantFeatureCounts.queue
  clusterOptions params.resources.quantFeatureCounts.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_FeatureCounts_quant", mode: 'symlink'
  input:
  tuple(val(flag), path(bam), path(bai))
  
  output:
  tuple(val(flag), path('*_featureCounts.txt'), path('*_featureCounts.txt.summary'), path('*.err'))

  shell:
  '''
  outname=!{flag}_featureCounts
  outfile=$outname'.txt'
  errorfile=$outname'.err'
  
  featureCounts -T !{params.resources.quantFeatureCounts.cpus} -p \
     -s !{params.quantFeatureCounts.strandedness} !{params.quantFeatureCounts.options}\
     -a !{params.quantFeatureCounts.annot} \
     -o $outfile !{bam} 2> $errorfile

  '''
}
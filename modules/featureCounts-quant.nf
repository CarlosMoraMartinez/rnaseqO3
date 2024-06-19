process quantFeatureCounts{
  label 'mg03_featureCounts_quant'
  conda params.quantFeatureCounts.conda
  cpus params.resources.quantFeatureCounts.cpus
  memory params.resources.quantFeatureCounts.mem
  queue params.resources.quantFeatureCounts.queue
  clusterOptions params.resources.quantFeatureCounts.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_featureCounts_quant", mode: 'symlink'
  input:
  tuple(val(flag), path(bam), path(bai))
  
  output:
  tuple(val(flag), path('*_featureCounts.txt'), path('*_featureCounts.txt.summary'), path('*_featureCounts.err'),
             path('*_featureCounts_byExon.txt'), path('*_featureCounts_byExon.txt.summary'), path('*_byExon.err'))

  shell:
  '''
  # Count by gene
  outname=!{flag}_featureCounts
  outfile=$outname'.txt'
  errorfile=$outname'.err'
  
  featureCounts -T !{params.resources.quantFeatureCounts.cpus} -p \
     -s !{params.quantFeatureCounts.strandedness} !{params.quantFeatureCounts.options}\
     -a !{params.quantFeatureCounts.annot} \
     -o $outfile !{bam} 2> $errorfile

  outname=!{flag}_featureCounts_byExon
  outfile=$outname'.txt'
  errorfile=$outname'.err'

  # Count by exon
  featureCounts -T !{params.resources.quantFeatureCounts.cpus} -p -f -t exon \
     -s !{params.quantFeatureCounts.strandedness} !{params.quantFeatureCounts.options}\
     -a !{params.quantFeatureCounts.annot} \
     -o $outfile !{bam} 2> $errorfile

  '''
}
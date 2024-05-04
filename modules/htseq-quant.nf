
process quantHTSeq{
  label 'mg03_htseq_quant'
  conda params.quantHTSeq.conda
  cpus params.resources.quantHTSeq.cpus
  memory params.resources.quantHTSeq.mem
  queue params.resources.quantHTSeq.queue
  clusterOptions params.resources.quantHTSeq.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_htseq_quant", mode: 'symlink'
  input:
  tuple(val(flag), path(bam), path(bai))
  
  output:
  tuple(val(flag), path('*_union.txt'), path('*_union.err'))

  shell:
  '''
  outname=!{flag}_htseq

  # Do all three different modes od HTSeq
  outfile_union=$outname'_union.txt'
  errorfile_union=$outname'_union.err'
  
  htseq-count -f bam -r pos \
    -m union !{params.quantHTSeq.options} \
    -a !{params.quantHTSeq.minaqual} \
    --stranded !{params.quantHTSeq.stranded} \
    --secondary-alignments !{params.quantHTSeq.secondary_alignments} \
    --supplementary-alignments !{params.quantHTSeq.supplementary_alignments} \
    !{bam} !{params.quantHTSeq.annot}> $outfile_union 2>$errorfile_union

  '''
}


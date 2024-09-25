
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
  tuple(val(mode), val(nonunique), val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val(flag), path('*.txt'), path('*.err'))

  shell:
  '''

  outname=!{flag}_!{sample_id}_htseq_!{mode}_!{nonunique}

  # Do all three different modes od HTSeq
  outfile_union=$outname'.txt'
  errorfile_union=$outname'.err'
  
  htseq-count -f bam -r pos --with-header -n !{params.resources.quantHTSeq.cpus} \
    -m !{mode} --nonunique !{nonunique} !{params.quantHTSeq.options} \
    -a !{params.quantHTSeq.minaqual} \
    --stranded !{params.quantHTSeq.stranded} \
    --secondary-alignments !{params.quantHTSeq.secondary_alignments} \
    --supplementary-alignments !{params.quantHTSeq.supplementary_alignments} \
    !{bam} !{params.quantHTSeq.annot}> $outfile_union 2>$errorfile_union
  
  # Clean sample names in the header
  sed -i !{params.quantHTSeq.cleaning_names_regex} *.txt
  '''
}


process quantStringtie2{
  label 'mg03_stringtie2_quant'
  conda params.quantStringtie2.conda
  cpus params.resources.quantStringtie2.cpus
  memory params.resources.quantStringtie2.mem
  queue params.resources.quantStringtie2.queue
  clusterOptions params.resources.quantStringtie2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_stringtie2_quant", mode: 'symlink'
  input:
  tuple(val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val(flag), val(sample_id), path('*_Stringtie'), path('*.err'))

  shell:
  '''
  outdir=!{sample_id}_!{flag}_Stringtie
  outgtf=!{sample_id}'.gtf'
  outtsv=!{sample_id}'.tsv'
  errorfile=$outdir'.err'
  
  stringtie !{bam} -p !{params.resources.quantStringtie2.cpus} -e !{params.quantStringtie2.options} \
            -G !{params.quantStringtie2.annot} \
            -o $outgtf \
            -A $outtsv \
            -b $outdir 2> $errorfile
  mv $outgtf $outtsv $outdir
  '''
}
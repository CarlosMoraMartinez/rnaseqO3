process doBBduk{
  label 'mg01_bbduk'
  conda params.doBBduk.conda
  cpus params.resources.doBBduk.cpus
  memory params.resources.doBBduk.mem
  queue params.resources.doBBduk.queue
  clusterOptions params.resources.doBBduk.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg01_bbduk", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(fastq))
  
  output:
  tuple(val(illumina_id), path('*.bbduk.fastq.gz'), path('*_bbduk.log'), path('*_bbduk.err'))

  shell:
  '''
  logfile=!{illumina_id}_bbduk.log
  errfile=!{illumina_id}_bbduk.err

  bbduk.sh in1=!{fastq[0]} \
           in2=!{fastq[1]}  \
           out1=!{illumina_id}_1.bbduk.fastq.gz  \
           out2=!{illumina_id}_2.bbduk.fastq.gz  \
           ref=!{params.doBBduk.illuminaclip} \
           threads=!{params.resources.doBBduk.cpus} \
           trimpolya=!{params.doBBduk.trimpolya}  \
           ktrim=r  \
           k=!{params.doBBduk.k} \
           mink=!{params.doBBduk.mink}  \
           hdist=!{params.doBBduk.hdist}  \
           tpe=f tbo=f !{params.doBBduk.extra_args} >$logfile 2>$errfile
  '''
}
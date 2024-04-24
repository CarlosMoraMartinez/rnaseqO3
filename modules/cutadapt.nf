process doCutadapt{
  label 'mg01_cutadapt'
  conda params.doCutadapt.conda
  cpus params.resources.doCutadapt.cpus
  memory params.resources.doCutadapt.mem
  queue params.resources.doCutadapt.queue
  clusterOptions params.resources.doCutadapt.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg01_cutadapt", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(fastq))

  output:
  tuple(val(illumina_id), path('*.cutadapt.fastq.gz'), path('*_cutadapt.log'), path('*_cutadapt.err'))

  shell:
  '''
  #Get adapter sequences from Fasta with '>Forward' and '>Reverse' entries
  FWFILT=$(grep -A 1 Forward !{params.doCutadapt.illuminaclip} | grep -v ">")
  RVFILT=$(grep -A 1 Reverse !{params.doCutadapt.illuminaclip} | grep -v ">")

  logfile=!{illumina_id}_cutadapt.log
  errfile=!{illumina_id}_cutadapt.err

  cutadapt \
    -j !{params.resources.doCutadapt.cpus} \
    -a $FWFILT \
    -A $RVFILT \
    --minimum-length !{params.doCutadapt.minlength} \
    -q !{params.doCutadapt.quality_cutoff_fw} \
    -Q !{params.doCutadapt.quality_cutoff_rv} \
    !{params.doCutadapt.extra_args} \
    -o !{illumina_id}_1.cutadapt.fastq.gz \
    -p !{illumina_id}_2.cutadapt.fastq.gz \
    !{fastq[0]} !{fastq[1]} >$logfile 2>$errfile
  '''
}
process picardMarkDuplicates{
  label 'mg04_picard_markduplicates'
  conda params.picardMarkDuplicates.conda
  cpus params.resources.picardMarkDuplicates.cpus
  memory params.resources.picardMarkDuplicates.mem
  queue params.resources.picardMarkDuplicates.queue
  clusterOptions params.resources.picardMarkDuplicates.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg04_picard_markduplicates", mode: 'symlink'
  input:
  tuple(val(flag), val(sample_id), path(bam), path(bai))
  
  output:
  tuple(val(flag), val(sample_id), path("*.markdups.bam"), path("*.markdups.bam.bai"), path("*.txt"))

  shell:
  '''
  newinputbam=!{sample_id}_!{flag}'.sorted.bam'
  newinputbai=$newinputbam.bai

  outreport=!{sample_id}-!{flag}_picard_markduplicates.txt
  outbam=!{sample_id}-!{flag}'.markdups.bam'
  errorfile=!{sample_id}-!{flag}_picard_markduplicates.err
  outbai=!{sample_id}-!{flag}'.markdups.bai'
  mem=$(echo !{params.resources.picardMarkDuplicates.mem} | sed "s/ GB/G/")

  mv !{bam} $newinputbam
  mv !{bai} $newinputbai

  java -Xms$mem -Xmx$mem \
              -jar !{params.software.picard_path} MarkDuplicates \
              --INPUT $newinputbam --OUTPUT $outbam !{params.picardMarkDuplicates.options} \
              --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
              --METRICS_FILE $outreport 2> $errorfile
  
  mv $outbai $outbam'.bai'
  '''
}
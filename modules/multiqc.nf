process multiQC{
  label 'mg11_multiqc'
  conda params.multiQC.conda
  cpus params.resources.multiQC.cpus
  memory params.resources.multiQC.mem
  queue params.resources.multiQC.queue
  clusterOptions params.resources.multiQC.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg11_multiqc", mode: 'copy'
  input:
    path yaml
    path all_reports

  output:
  path("multiqc_report.html")

  shell:
  '''
  mv !{yaml} multiqc_config.yaml
  multiqc --filename multiqc_report.html .  
  '''
}
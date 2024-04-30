process mergeStringtie2{
  label 'mg03_stringtie2_merge'
  conda params.mergeStringtie2.conda
  cpus params.resources.mergeStringtie2.cpus
  memory params.resources.mergeStringtie2.mem
  queue params.resources.mergeStringtie2.queue
  clusterOptions params.resources.mergeStringtie2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_stringtie2_merge", mode: 'symlink'
  input:
  tuple(val(flag), path(all_sample_counts))
  
  output:
  tuple(val(flag), path('*gene_count_matrix.csv'), path('*transcript_count_matrix.csv'), path('*.err'))

  shell:
  '''
  out_genecounts=!{flag}_Stringtie_gene_count_matrix.csv
  out_transccounts=!{flag}_Stringtie_transcript_count_matrix.csv
  errorfile=!{flag}_Stringtie_merge.err
  
  prepDE.py -g $out_genecounts -t $out_transccounts \
            -l !{params.resources.mergeStringtie2.readlength} 2>$errorfile
  
  sed -i "s/_!{flag}_Stringtie//g" $out_genecounts 
  sed -i "s/_!{flag}_Stringtie//g" $out_transccounts 

  '''
}
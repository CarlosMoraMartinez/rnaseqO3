process buildindexSalmon{
  label 'mg00_salmon_buildindex'
  conda params.buildindexSalmon.conda
  cpus params.resources.buildindexSalmon.cpus
  memory params.resources.buildindexSalmon.mem
  queue params.resources.buildindexSalmon.queue
  clusterOptions params.resources.buildindexSalmon.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_salmon_buildindex", mode: 'symlink'
  input:
  val(index_dir)
  val(mode)
  path(genome_fasta)
  path(transcripts_fasta)
  path(annot)
  
  output:
  path("${index_dir}_${mode}")

  shell:
  if(mode == "transcripts_only")
  '''
  outdir=!{index_dir}'_'!{mode}

  salmon index -t !{transcripts_fasta} -i $outdir -k !{params.buildindexSalmon.k} !{params.buildindexSalmon.options}
  '''

  else if(mode == "genome_decoy")
  '''
  outdir=!{index_dir}'_'!{mode}
  
  cat !{transcripts_fasta} !{genome_fasta} > transcripts_genomedecoy.fasta
  grep ">" !{genome_fasta} | sed "s/>//" | cut -f1 -d' ' > decoys.txt

  salmon index -t transcripts_genomedecoy.fasta -i $outdir --decoys decoys.txt -k !{params.buildindexSalmon.k} !{params.buildindexSalmon.options}

  rm transcripts_genomedecoy.fasta
 
  '''
}
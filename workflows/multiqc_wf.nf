include { multiQC } from '../modules/multiqc'


workflow MULTIQC{
  take:
  ch_fastqc
  ch_picard_rnametrics
  ch_alignment_markdups
  ch_rnaseqc
  ch_hisat2_result
  ch_star_result
  ch_star_2ndpass_result
  ch_salmon_aln_result
  ch_salmon_merged
  ch_kallisto_result
  ch_fcounts_results
  ch_htseq_results
  ch_subread_result
  ch_bbmap_result
  ch_fastq_processed_paired

  main:

  //prepare MultiQC process input
  all_reports = ch_fastqc.collect().ifEmpty([])
  .view()
  //trim_qc = ch_fastq_processed.map{it -> it[4]}.collect().ifEmpty([])
  //bowtie2_err = ch_alignment_output.map{it -> it[2]}.collect().ifEmpty([])
  //kraken_err = ch_kraken2_output.map{it -> it[2]}.collect().ifEmpty([])
  //bracken_err = ch_bracken_output.map{it -> it[4]}.collect().ifEmpty([])
  //metaphlan_err = ch_metaphlan_output.map{it -> it[1]}.collect().ifEmpty([])
  //megahit_err = ch_megahit_output.map{it -> it[4]}.collect().ifEmpty([])
  //ch_metaquast_tsv = ch_metaquast_output.map{it -> it[3]}.collect().ifEmpty([])

  //Call multiQC process
  multiQC(params.multiQC.configyaml,
            all_reports
            )
  ch_multiqc_out = multiQC.out //.map{it -> it[1]}
emit:
ch_multiqc_out
}
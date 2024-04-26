


include { alignHISAT2 } from '../modules/hisat2-align'
include { alignSTAR } from '../modules/star-align'

workflow ALIGN {
  take: ch_fastq_processed_paired
  main:
 //Align using HISAT2
  ch_hisat2_result = Channel.from([])
  ch_hisat2_bam = Channel.from([])
  if(params.alignHISAT2.do_hisat2){ 
    alignHISAT2(ch_fastq_processed_paired)
     ch_hisat2_result = alignHISAT2.out
     //.view{ "HISAT2 full result: $it" }
     ch_hisat2_bam = ch_hisat2_result.map{it -> tuple("HISAT2", it[0], it[1], it[2])}
     .view{ "HISAT2 BAM only: $it" }
  }

  ch_star_result = Channel.from([])
  ch_star_bam = Channel.from([])
  if(params.alignSTAR.do_star){ 
    alignSTAR(ch_fastq_processed_paired)
     ch_star_result = alignSTAR.out
     //.view{ "HISAT2 full result: $it" }
     ch_star_bam = ch_star_result.map{it -> tuple("STAR", it[0], it[1])}
     .view{ "STAR BAM only: $it" }
  }

  emit:
  ch_hisat2_result
  ch_hisat2_bam
  ch_star_result
  ch_star_bam
}
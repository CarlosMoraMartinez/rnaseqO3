


include { buildindexHISAT2 } from '../modules/hisat2-buildindex'
include { alignHISAT2 } from '../modules/hisat2-align'
include { alignSTAR } from '../modules/star-align'
include { buildindexSTAR } from '../modules/star-buildindex'
include { alignSTAR2ndPass } from '../modules/star-align-2ndpass'

workflow ALIGN {
  take: ch_fastq_processed_paired
  main:
  //Align using HISAT2
  ch_hisat2_result = Channel.from([])
  ch_hisat2_bam = Channel.from([])

  if(params.alignHISAT2.do_hisat2){ 
    if(params.buildindexHISAT2.do_index){
      buildindexHISAT2(
        params.buildindexHISAT2.index_dir,
        params.buildindexHISAT2.index_name,
        params.buildindexHISAT2.fasta,
        params.buildindexHISAT2.annot)
      ch_hisat2_index = buildindexHISAT2.out
        .map{it -> it[0]}
        .view{ "HISAT2 index created: $it" }
    }else{
      ch_hisat2_index = params.alignHISAT2.index
    }
    alignHISAT2(ch_hisat2_index, ch_fastq_processed_paired)
    ch_hisat2_result = alignHISAT2.out
     //.view{ "HISAT2 full result: $it" }
    ch_hisat2_bam = ch_hisat2_result.map{it -> tuple("HISAT2", it[0], it[1], it[2])}
      .view{ "HISAT2 BAM only: $it" }
  }

   //Align using STAR
  ch_star_result = Channel.from([])
  ch_star_bam = Channel.from([])
  ch_star_2ndpass_result = Channel.from([])
  ch_star_2ndpass_bam = Channel.from([])

  if(params.alignSTAR.do_star){ 
    if(params.buildindexSTAR.do_index){
      buildindexSTAR(params.buildindexSTAR.index_name,
        params.buildindexSTAR.fasta,
        params.buildindexSTAR.annot)
      ch_star_index = buildindexSTAR.out
       //.view{ "STAR index created: $it" }
    }else{
      ch_star_index = params.alignSTAR.index
    }
    alignSTAR(ch_star_index, ch_fastq_processed_paired)
     ch_star_result = alignSTAR.out
     //.view{ "STAR full result: $it" }
     ch_star_bam = ch_star_result.map{it -> tuple("STAR", it[0], it[1], it[2])}
     .view{ "STAR BAM only: $it" }

  if(params.alignSTAR2ndPass.do_star){
    ch_star2ndpass_input_SJ = ch_star_result.map{it -> tuple(it[3])}
       .collect()
       .view{ "STAR splice junctions collect: $it" }
    alignSTAR2ndPass(ch_star_index, ch_fastq_processed_paired, ch_star2ndpass_input_SJ)
    ch_star_2ndpass_result = alignSTAR2ndPass.out
    //.view{ "STAR 2nd PASS full result: $it" }
     ch_star_2ndpass_bam = ch_star_2ndpass_result.map{it -> tuple("STAR2", it[0], it[1], it[2])}
    .view{ "STAR 2nd PASS BAM only: $it" }
  }

    
  }

  emit:
  ch_hisat2_result
  ch_hisat2_bam
  ch_star_result
  ch_star_bam
  ch_star_2ndpass_result
  ch_star_2ndpass_bam
}
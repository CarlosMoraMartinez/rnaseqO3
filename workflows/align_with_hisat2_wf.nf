include { buildindexHISAT2 } from '../modules/hisat2-buildindex'
include { alignHISAT2 } from '../modules/hisat2-align'

workflow ALIGN_WITH_HISAT2 {
  take: ch_fastq_processed_paired

  main:

  if(params.buildindexHISAT2.do_index){
    buildindexHISAT2(
      params.buildindexHISAT2.index_dir,
      params.buildindexHISAT2.index_name,
      params.buildindexHISAT2.fasta,
      params.buildindexHISAT2.annot)
    ch_hisat2_index = buildindexHISAT2.out
      .map{it -> it[0]}
      //.view{ "HISAT2 index created: $it" }
  }else{
    ch_hisat2_index = params.alignHISAT2.index
  }
  alignHISAT2(ch_hisat2_index, ch_fastq_processed_paired)
  ch_hisat2_result = alignHISAT2.out
   //.view{ "HISAT2 full result: $it" }
  ch_hisat2_bam = ch_hisat2_result.map{it -> tuple("HISAT2", it[0], it[1], it[2])}
    //.view{ "HISAT2 BAM only: $it" }

  emit:
  ch_hisat2_result
  ch_hisat2_bam
}
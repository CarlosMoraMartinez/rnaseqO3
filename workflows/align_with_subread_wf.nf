include { buildindexSubread } from '../modules/subread-buildindex'
include { alignSubread } from '../modules/subread-align'

workflow ALIGN_WITH_SUBREAD {
  take: ch_fastq_processed_paired

  main:

  if(params.buildindexSubread.do_index){
    buildindexSubread(
      params.buildindexSubread.index_dir,
      params.buildindexSubread.index_name,
      params.buildindexSubread.fasta
      )
    ch_Subread_index = buildindexSubread.out
      .map{it -> it[0]}
      //.view{ "Subread index created: $it" }
  }else{
    ch_Subread_index = params.alignSubread.index_prefix
  }

  alignSubread(ch_Subread_index, ch_fastq_processed_paired)
  ch_subread_result = alignSubread.out
   //.view{ "Subread full result: $it" }
  ch_subread_bam = ch_subread_result.map{it -> tuple("Subread", it[0], it[1], it[2])}
    //.view{ "Subread BAM only: $it" }

  emit:
  ch_subread_result
  ch_subread_bam
}
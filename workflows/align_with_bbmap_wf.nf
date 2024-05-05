include { buildindexBBMap } from '../modules/bbmap-buildindex'
include { alignBBMap } from '../modules/bbmap-align'

workflow ALIGN_WITH_BBMAP {
  take: ch_fastq_processed_paired

  main:

  if(params.buildindexBBMap.do_index){
    buildindexBBMap(params.buildindexBBMap.fasta)
    ch_BBMap_index = buildindexBBMap.out
  }else{
    ch_BBMap_index = params.alignBBMap.index
  }

  alignBBMap(ch_BBMap_index, ch_fastq_processed_paired)

  ch_bbmap_result = alignBBMap.out
   .view{ "BBMap full result: $it" }
  ch_bbmap_bam = ch_bbmap_result.map{it -> tuple("BBMap", it[0], it[1], it[2])}
    .view{ "BBMap BAM only: $it" }

  emit:
  ch_bbmap_result
  ch_bbmap_bam
}
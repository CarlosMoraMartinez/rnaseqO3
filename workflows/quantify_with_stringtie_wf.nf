include { quantStringtie2 } from '../modules/stringtie2-quant'
include { mergeStringtie2 } from '../modules/stringtie2-merge'

workflow QUANTIFY_WITH_STRINGTIE {
  take: 
  ch_alignment_all

  main:

  // Call Stringtie2
  quantStringtie2(ch_alignment_all)
  ch_stringtie_results = quantStringtie2.out
  //.view{ "Stringtie results: $it" }

  // Group results by aligner
  ch_stringtie_results_grouped = ch_stringtie_results
    .map{it -> tuple(it[0], tuple(it[2]))}
    .groupTuple(by: 0)
    .map{it -> tuple(it[0], it[1].flatten())}
    //.view{ "Stringtie results grouped: $it" }

  // Merge results from all samples, by aligner
  if(params.mergeStringtie2.do_merge){
  mergeStringtie2(ch_stringtie_results_grouped)
  ch_stringtie_results_merged = mergeStringtie2.out
    //.view{ "Stringtie merge prepDE.py results: $it" }
  } else {
    ch_stringtie_results_merged = Channel.from([])
  }
  emit:
  ch_stringtie_results
  ch_stringtie_results_merged
}
include { ALIGN_WITH_HISAT2 } from './align_with_hisat2_wf'
include { ALIGN_WITH_STAR } from './align_with_star_wf'

include { quantStringtie2 } from '../modules/stringtie2-quant'
include { mergeStringtie2 } from '../modules/stringtie2-merge'
include { quantSalmon } from '../modules/salmon-quant'
include { quantBamSalmon } from '../modules/salmon-quant-frombam'
include { mergeSalmon } from '../modules/salmon-merge'
include { mapdecoySalmontools } from '../modules/salmontools-mapdecoy'
include { buildindexSalmon } from '../modules/salmon-buildindex'

workflow QUANTIFY_WITH_STRINGTIE {
  take: 
  ch_hisat2_bam
  ch_star_bam
  ch_star_2ndpass_bam

  main:

  //First concatenate all alignment channels
  ch_alignment_all = ch_hisat2_bam
    .concat(ch_star_bam)
    .concat(ch_star_2ndpass_bam)
    //.view{ "All alignments concat: $it" }

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
  mergeStringtie2(ch_stringtie_results_grouped)
  ch_stringtie_results_merged = mergeStringtie2.out
    //.view{ "Stringtie merge prepDE.py results: $it" }

  emit:
  ch_stringtie_results
  ch_stringtie_results_merged
}
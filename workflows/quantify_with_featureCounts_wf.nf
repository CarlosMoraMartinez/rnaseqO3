include { quantFeatureCounts } from '../modules/featureCounts-quant'

workflow QUANTIFY_WITH_FEATURECOUNTS {
  take: 
  ch_alignment_all

  main:
  
  ch_alignment_all_grouped = ch_alignment_all
    .groupTuple(by: 0)
    .map{it -> tuple(it[0], it[2].flatten(), it[3].flatten())}
    //.view{ "featureCounts input grouped: $it" }

  // Call FeatureCounts only once with all samples
  quantFeatureCounts(ch_alignment_all_grouped)
  ch_fcounts_results = quantFeatureCounts.out
    //.view{ "featureCounts results: $it" }


  emit:
  ch_fcounts_results
}
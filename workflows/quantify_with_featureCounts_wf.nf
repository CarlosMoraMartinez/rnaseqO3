include { quantFeatureCounts } from '../modules/featureCounts-quant'

workflow QUANTIFY_WITH_FEATURECOUNTS {
  take: 
  ch_hisat2_bam
  ch_star_bam
  ch_star_2ndpass_bam
  ch_subread_bam
  ch_bbmap_bam

  main:

  //First concatenate all alignment channels
  ch_alignment_all = ch_hisat2_bam
    .concat(ch_subread_bam)
    .concat(ch_bbmap_bam)
    //.concat(ch_star_bam)
    //.concat(ch_star_2ndpass_bam)
    //.view{ "All alignments concat: $it" }
  
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
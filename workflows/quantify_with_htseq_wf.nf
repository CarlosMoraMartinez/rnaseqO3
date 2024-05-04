include { quantHTSeq } from '../modules/htseq-quant'

workflow QUANTIFY_WITH_HTSEQ {
  take: 
  ch_hisat2_bam
  ch_star_bam
  ch_star_2ndpass_bam
  ch_subread_bam

  main:

  //First concatenate all alignment channels
  ch_alignment_all = ch_hisat2_bam
    .concat(ch_subread_bam)
    .concat(ch_star_bam)
    .concat(ch_star_2ndpass_bam)
    //.view{ "All alignments concat: $it" }
  
  ch_alignment_all_grouped = ch_alignment_all
    .groupTuple(by: 0)
    .map{it -> tuple(it[0], it[2].flatten(), it[3].flatten())}
    //.view{ "htseq input grouped: $it" }

  // Run HT-seq in many different modes

 
  ch_htseq_mode = Channel.fromList(params.quantHTSeq.mode_list.tokenize(','))
                  .combine(Channel.fromList(params.quantHTSeq.nonunique_list.tokenize(',')))
                  .view{ "htseq modes combinations: $it" }

  ch_alignment_all_grouped2 = ch_htseq_mode.combine(ch_alignment_all_grouped)
    //.view{ "htseq input COMBINED: $it" }

  // Call HTSeq only once with all samples
  quantHTSeq(ch_alignment_all_grouped2)
  ch_htseq_results = quantHTSeq.out
    //.view{ "htseq results: $it" }


  emit:
  ch_htseq_results
}
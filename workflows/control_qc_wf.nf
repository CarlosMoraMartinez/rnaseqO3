include { buildRefflatFromGTF } from '../modules/build_refflat_from_gtf'
include { buildIntervalListFromBed } from '../modules/picard_intervallist_from_bed'
include { picardMarkDuplicates } from '../modules/picard_mark_duplicates'
include { picardRNASeqMetrics } from '../modules/picard_rnaseqmetrics'

workflow CONTROL_QC{
  take:
  ch_hisat2_result
  ch_hisat2_bam
  ch_star_result
  ch_star_bam
  ch_star_2ndpass_result
  ch_star_2ndpass_bam
  ch_subread_result
  ch_subread_bam
  ch_bbmap_result
  ch_bbmap_bam
  ch_alignment_all
  ch_stringtie_results
  ch_stringtie_results_merged
  ch_salmon_result
  ch_salmon_aln_result
  ch_salmon_merged
  ch_kallisto_result
  ch_fcounts_results
  ch_htseq_results

  main:

  // 0. Picard markDuplicates
  ch_alignment_2_qual= ch_htseq_results
  if(params.picardMarkDuplicates.do){
    picardMarkDuplicates(ch_alignment_all)
    ch_alignment_markdups = picardMarkDuplicates.out
    ch_alignment_2_qual = ch_alignment_markdups
    .map{it -> tuple(it[0], it[1], it[2], it[3])}
    //.view{ "MarkDuplicates output: $it" }
  }
  
  // 1. Picard RNASeq metrics
  ch_refflat = params.picardRNASeqMetrics.refflat
  ch_rrna_intervallist = params.picardRNASeqMetrics.rRNA_interval_list
  ch_picard_rnametrics = Channel.from([])
  if(params.picardRNASeqMetrics.do){

    //Build refflat file from GTF annotation if necessary
    if(params.buildRefflatFromGTF.do){
      buildRefflatFromGTF(params.buildRefflatFromGTF.genome_fasta, 
                          params.buildRefflatFromGTF.annot)
      ch_refflat = buildRefflatFromGTF.out
      .view{ "REFFLAT file built: $it" }
    }

    //Build interval_list file from BED with rRNA regions if necessary
    if(params.buildIntervalListFromBed.do){
      buildIntervalListFromBed(params.buildIntervalListFromBed.genome_fasta,
                              params.buildIntervalListFromBed.bed)
      ch_rrna_intervallist = buildIntervalListFromBed.out
      .view{ "INTERVAL_LIST file built: $it" }
    }

    picardRNASeqMetrics(ch_refflat, ch_rrna_intervallist, ch_alignment_2_qual)
    ch_picard_rnametrics = picardRNASeqMetrics.out
    //.view{ "Picard RNASeqMetrics result: $it" }
  }

  emit:
  ch_alignment_markdups
  ch_picard_rnametrics
}
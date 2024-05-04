include { ALIGN_WITH_HISAT2 } from './align_with_hisat2_wf'
include { ALIGN_WITH_SUBREAD } from './align_with_subread_wf'
include { ALIGN_WITH_BBMAP } from './align_with_bbmap_wf'
include { ALIGN_WITH_STAR } from './align_with_star_wf'
include { QUANTIFY_WITH_STRINGTIE } from './quantify_with_stringtie_wf'
include { QUANTIFY_WITH_FEATURECOUNTS } from './quantify_with_featureCounts_wf'
include { QUANTIFY_WITH_HTSEQ } from './quantify_with_htseq_wf'
include { QUANTIFY_WITH_SALMON } from './quantify_with_salmon_wf'
include { QUANTIFY_WITH_KALLISTO } from './quantify_with_kallisto_wf'

workflow ALIGN_ALL {
  take: ch_fastq_processed_paired
  main:

  //Align using HISAT2
  ch_hisat2_result = Channel.from([])
  ch_hisat2_bam = Channel.from([])

  if(params.workflows.do_hisat2){ 
    ALIGN_WITH_HISAT2(ch_fastq_processed_paired)
    ch_hisat2_result = ALIGN_WITH_HISAT2.out.ch_hisat2_result
    ch_hisat2_bam = ALIGN_WITH_HISAT2.out.ch_hisat2_bam
  } // end HISAT2

   //Align using Subread
  ch_subread_result = Channel.from([])
  ch_subread_bam = Channel.from([])

  if(params.workflows.do_subread){ 
    ALIGN_WITH_SUBREAD(ch_fastq_processed_paired)
    ch_subread_result = ALIGN_WITH_SUBREAD.out.ch_subread_result
    ch_subread_bam = ALIGN_WITH_SUBREAD.out.ch_subread_bam
  } // end Subread

  //Align using BBMap
  ch_bbmap_result = Channel.from([])
  ch_bbmap_bam = Channel.from([])

  if(params.workflows.do_bbmap){ 
    ALIGN_WITH_BBMAP(ch_fastq_processed_paired)
    ch_bbmap_result = ALIGN_WITH_BBMAP.out.ch_bbmap_result
    ch_bbmap_bam = ALIGN_WITH_BBMAP.out.ch_bbmap_bam
  } // end BBMap


  //Align using STAR
  ch_star_result = Channel.from([])
  ch_star_bam = Channel.from([])
  ch_star_bam_bytranscript = Channel.from([])
  ch_star_2ndpass_result = Channel.from([])
  ch_star_2ndpass_bam = Channel.from([])
  ch_star_2ndpass_bam_bytranscript = Channel.from([])

  if(params.workflows.do_star){ 

    ALIGN_WITH_STAR(ch_fastq_processed_paired)

    ch_star_result = ALIGN_WITH_STAR.out.ch_star_result
    ch_star_bam = ALIGN_WITH_STAR.out.ch_star_bam
    ch_star_bam_bytranscript = ALIGN_WITH_STAR.out.ch_star_bam_bytranscript
    ch_star_2ndpass_result = ALIGN_WITH_STAR.out.ch_star_2ndpass_result
    ch_star_2ndpass_bam = ALIGN_WITH_STAR.out.ch_star_2ndpass_bam
    ch_star_2ndpass_bam_bytranscript = ALIGN_WITH_STAR.out.ch_star_2ndpass_bam_bytranscript
  }// end STAR

  // Quantify using STRINGTIE
  ch_stringtie_results_merged = Channel.from([])
  ch_stringtie_results = Channel.from([])
  if(params.workflows.do_stringtie){

    QUANTIFY_WITH_STRINGTIE(
      ch_hisat2_bam, 
      ch_star_bam, 
      ch_star_2ndpass_bam,
      ch_subread_bam,
      ch_bbmap_bam
      
    )

    ch_stringtie_results_merged = QUANTIFY_WITH_STRINGTIE.out.ch_stringtie_results_merged
    ch_stringtie_results = QUANTIFY_WITH_STRINGTIE.out.ch_stringtie_results
  } // end STRINGTIE

  // Quantify using featureCounts
  ch_fcounts_results = Channel.from([])
  if(params.workflows.do_featureCounts){

    QUANTIFY_WITH_FEATURECOUNTS(
      ch_hisat2_bam, 
      ch_star_bam, 
      ch_star_2ndpass_bam,
      ch_subread_bam,
      ch_bbmap_bam
    )

    ch_fcounts_results = QUANTIFY_WITH_FEATURECOUNTS.out.ch_fcounts_results
  } // end featureCounts

  // Quantify using HTSeq
  ch_htseq_results = Channel.from([])
  if(params.workflows.do_htseq){

    QUANTIFY_WITH_HTSEQ(
      ch_hisat2_bam, 
      ch_star_bam, 
      ch_star_2ndpass_bam,
      ch_subread_bam,
      ch_bbmap_bam
    )

    ch_htseq_results = QUANTIFY_WITH_HTSEQ.out.ch_htseq_results
  } // end HTSeq

  // Quantify using   SALMON
  ch_salmon_result = Channel.from([])
  ch_salmon_aln_result = Channel.from([])
  ch_salmon_merged = Channel.from([])

  if(params.workflows.do_salmon){

    QUANTIFY_WITH_SALMON ( 
      ch_fastq_processed_paired,
      ch_star_bam_bytranscript,
      ch_star_2ndpass_bam_bytranscript
    )

    ch_salmon_result = QUANTIFY_WITH_SALMON.out.ch_salmon_result
    ch_salmon_aln_result = QUANTIFY_WITH_SALMON.out.ch_salmon_aln_result
    ch_salmon_merged = QUANTIFY_WITH_SALMON.out.ch_salmon_merged
  } // end SALMON

  // Quantify with Kallisto
  ch_kallisto_result = Channel.from([])
  if(params.workflows.do_kallisto){
    QUANTIFY_WITH_KALLISTO(ch_fastq_processed_paired)
    
    ch_kallisto_result = QUANTIFY_WITH_KALLISTO.out.ch_kallisto_result
  }

  emit:
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
  ch_stringtie_results
  ch_stringtie_results_merged
  ch_salmon_result
  ch_salmon_aln_result
  ch_salmon_merged
  ch_kallisto_result
  ch_fcounts_results
  ch_htseq_results
}
include { ALIGN_WITH_HISAT2 } from './align_with_hisat2_wf'
include { ALIGN_WITH_STAR } from './align_with_star_wf'
include { QUANTIFY_WITH_STRINGTIE } from './quantify_with_stringtie_wf'
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
      ch_star_2ndpass_bam
    )

    ch_stringtie_results_merged = QUANTIFY_WITH_STRINGTIE.out.ch_stringtie_results_merged
    ch_stringtie_results = QUANTIFY_WITH_STRINGTIE.out.ch_stringtie_results
  } // end STRINGTIE

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
  ch_stringtie_results
  ch_stringtie_results_merged
  ch_salmon_result
  ch_salmon_aln_result
  ch_salmon_merged
  ch_kallisto_result
}
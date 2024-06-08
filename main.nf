include { CLEANFASTQ } from './workflows/cleanfastq_wf.nf'
include { ALIGN_ALL } from './workflows/align_all_wf.nf'
include { CONTROL_QC } from './workflows/control_qc_wf.nf'
include { MULTIQC } from './workflows/multiqc_wf.nf'

workflow {

  
   ch_rawfastq = Channel.fromFilePairs(params.raw_fastq)
     .view{"FilePairs input: $it"}

   if(params.workflows.doCleanFastq){
      //Call clean fastq workflow
      CLEANFASTQ(ch_rawfastq)

      //Get outputs (not all of them are used)
      ch_fastq_processed  = CLEANFASTQ.out.ch_fastq_processed
      ch_fastq_processed_paired = CLEANFASTQ.out.ch_fastq_processed_paired
      ch_fastqc = CLEANFASTQ.out.ch_fastqc
   }else{
      println "Skipping CLEANFASTQ. Taking raw fastq as final fastq."
      ch_fastq_processed_paired = ch_rawfastq
      ch_fastq_processed  = Channel.from([])
      ch_fastqc = Channel.from([])
   }

   // Workflow to perform alignment and quantification with all programs (depending on config)
   ALIGN_ALL(ch_fastq_processed_paired)

   ch_hisat2_result = ALIGN_ALL.out.ch_hisat2_result
   ch_hisat2_bam  = ALIGN_ALL.out.ch_hisat2_bam
   ch_star_result = ALIGN_ALL.out.ch_star_result
   ch_star_bam = ALIGN_ALL.out.ch_star_bam
   ch_star_2ndpass_result = ALIGN_ALL.out.ch_star_2ndpass_result
   ch_star_2ndpass_bam = ALIGN_ALL.out.ch_star_2ndpass_bam
   ch_subread_result = ALIGN_ALL.out.ch_subread_result
   ch_subread_bam = ALIGN_ALL.out.ch_subread_bam
   ch_bbmap_result = ALIGN_ALL.out.ch_bbmap_result
   ch_bbmap_bam = ALIGN_ALL.out.ch_bbmap_bam
   ch_alignment_all = ALIGN_ALL.out.ch_alignment_all
   ch_stringtie_results = ALIGN_ALL.out.ch_stringtie_results
   ch_stringtie_results_merged  = ALIGN_ALL.out.ch_stringtie_results_merged
   ch_salmon_result = ALIGN_ALL.out.ch_salmon_result
   ch_salmon_aln_result = ALIGN_ALL.out.ch_salmon_aln_result
   ch_salmon_merged = ALIGN_ALL.out.ch_salmon_merged
   ch_kallisto_result = ALIGN_ALL.out.ch_kallisto_result
   ch_fcounts_results = ALIGN_ALL.out.ch_fcounts_results
   ch_htseq_results = ALIGN_ALL.out.ch_htseq_results

   // Workflow to calculate quality metrics
   CONTROL_QC(
      ch_hisat2_result,
      ch_hisat2_bam,
      ch_star_result,
      ch_star_bam,
      ch_star_2ndpass_result,
      ch_star_2ndpass_bam,
      ch_subread_result,
      ch_subread_bam,
      ch_bbmap_result,
      ch_bbmap_bam,
      ch_alignment_all,
      ch_stringtie_results,
      ch_stringtie_results_merged,
      ch_salmon_result,
      ch_salmon_aln_result,
      ch_salmon_merged,
      ch_kallisto_result,
      ch_fcounts_results,
      ch_htseq_results
   )
   ch_alignment_markdups = CONTROL_QC.out.ch_alignment_markdups
   ch_picard_rnametrics = CONTROL_QC.out.ch_picard_rnametrics
   ch_rnaseqc = CONTROL_QC.out.ch_rnaseqc


   ////Call MultiQC workflow

   if(params.workflows.doMultiQC){
      MULTIQC(
        ch_fastqc,
        ch_picard_rnametrics,
        ch_alignment_markdups,
        ch_rnaseqc,
        ch_hisat2_result,
        ch_star_result,
        ch_star_2ndpass_result,
        ch_salmon_aln_result,
        ch_salmon_merged,
        ch_kallisto_result,
        ch_fcounts_results,
        ch_htseq_results,
        ch_subread_result,
        ch_bbmap_result,
        ch_fastq_processed_paired
      )
      ch_multiqc_out = MULTIQC.out.ch_multiqc_out
   }else{
      print "Skipping MULTIQC."
      ch_multiqc_out = Channel.from([])
   }
}
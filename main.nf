

include { CLEANFASTQ } from './workflows/cleanfastqwf.nf'
include { ALIGN } from './workflows/alignwf.nf'

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
      print "Skipping CLEANFASTQ. Taking raw fastq as final fastq."
      ch_fastq_processed_paired = ch_rawfastq
      ch_fastq_processed  = Channel.from([])
      ch_fastqc = Channel.from([])
    }

   ALIGN(ch_fastq_processed_paired)

   ////Call kraken workflow
   //if(params.workflows.doKraken2Bracken){
   //   KRAKEN2BRACKEN(ch_fastq_filtered)
   //   //Get outputs (not all of them are used)
   //   ch_kraken2_output = KRAKEN2BRACKEN.out.ch_kraken2_output
   //   ch_bracken_output = KRAKEN2BRACKEN.out.ch_bracken_output
   //   ch_transform2mpa_output = KRAKEN2BRACKEN.out.ch_transform2mpa_output
   //   ch_combineMpa_output = KRAKEN2BRACKEN.out.ch_combineMpa_output
   //   ch_krona_output = KRAKEN2BRACKEN.out.ch_krona_output
   //}else{
   //   print "Skipping KRAKEN2BRACKEN."
   //   ch_kraken2_output = Channel.from([])
   //   ch_bracken_output = Channel.from([])
   //   ch_transform2mpa_output = Channel.from([])
   //   ch_combineMpa_output = Channel.from([])
   //   ch_krona_output = Channel.from([])
   //}
//
  ////Call Humann3 workflow
   //if(params.workflows.doHumann3){
   //   HUMANN3(ch_fastq_filtered)
   //   ch_metaphlan = HUMANN3.out.ch_metaphlan
   //   ch_metaphlan_merged = HUMANN3.out.ch_metaphlan
   //   ch_humann3 = HUMANN3.out.ch_humann3
   //}else{
   //   ch_metaphlan = Channel.from([])
   //   ch_metaphlan_merged = Channel.from([])
   //   ch_humann3 = Channel.from([])
   //}
//
   ////Call Assembly workflow
   //if(params.workflows.doAssembly){
   //   ASSEMBLY(ch_fastq_filtered)
   //   ch_spades_output = ASSEMBLY.out.ch_spades_output
   //   ch_megahit_output = ASSEMBLY.out.ch_megahit_output
   //   ch_metaquast_output = ASSEMBLY.out.ch_metaquast_output
   //}else{
   //   ch_spades_output = Channel.from([])
   //   ch_megahit_output = Channel.from([])
   //   ch_metaquast_output = Channel.from([])
   //}
   ////Call MultiQC workflow
//
   //if(params.workflows.doMultiQC){
   //   MULTIQC(
   //     ch_fastqc,
   //     ch_fastq_processed,
   //     ch_alignment_output,
   //     ch_kraken2_output,
   //     ch_bracken_output,
   //     ch_metaphlan,
   //     ch_megahit_output,
   //     ch_metaquast_output
   //   )
   //   ch_multiqc_out = MULTIQC.out.ch_multiqc_out
   //}else{
   //   print "Skipping MULTIQC."
   //   ch_multiqc_out = Channel.from([])
   //}
}
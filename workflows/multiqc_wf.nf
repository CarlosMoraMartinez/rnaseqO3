include { multiQC } from '../modules/multiqc'


workflow MULTIQC{
  take:
  ch_fastqc
  ch_picard_rnametrics
  ch_alignment_markdups
  ch_rnaseqc
  ch_hisat2_result
  ch_star_result
  ch_star_2ndpass_result
  ch_salmon_aln_result
  ch_salmon_result
  ch_kallisto_result
  ch_fcounts_results
  ch_htseq_results
  ch_subread_result
  ch_bbmap_result
  ch_fastq_processed

  main:

  //prepare MultiQC process input

  //FastQC 
  all_reports = ch_fastqc.collect().ifEmpty([])

  // Trimming log
  if(params.workflows.do_trim && params.workflows.trimming_tool == "trimmomatic" ){ //With trimmomatic
     all_reports = ch_fastq_processed.map{it -> it[4]}.collect().ifEmpty([])
     .concat(all_reports).flatten().collect()   
  }else if (params.workflows.do_trim && params.workflows.trimming_tool == "cutadapt" ){ //With cutadapt
     all_reports = ch_fastq_processed.map{it -> it[2]}.collect().ifEmpty([])
     .concat(all_reports).flatten().collect()
  }else if (params.workflows.do_trim && params.workflows.trimming_tool == "bbduk" ){ //With bbduk
     all_reports = ch_fastq_processed.map{it -> it[3]}.collect().ifEmpty([])
     .concat(all_reports).flatten().collect()
  }
  //all_reports.view{"ch_fastq_processed: $it"}

  // HISAT2
  all_reports = ch_hisat2_result.map{it -> it[3..4]}.collect().ifEmpty([])
  //.view{"ch_fastq_processed: $it"}
  .concat(all_reports).flatten().collect()

  // STAR
  all_reports = ch_star_result.map{it -> it[6, 9..11]}.collect().ifEmpty([])
  //.view{"star: $it"}
  .concat(all_reports).flatten().collect()

  all_reports = ch_star_2ndpass_result.map{it -> it[6, 9..11]}.collect().ifEmpty([])
  //.view{"star 2: $it"}
  .concat(all_reports).flatten().collect()

  // Subread
  all_reports = ch_subread_result.map{it -> it[4, 5]}.ifEmpty([])
  //.view{"Subread: $it"}
  .concat(all_reports).flatten().collect()

  // BBMap
  all_reports = ch_bbmap_result.map{it -> it[3..5]}.flatten().collect().ifEmpty([])
  //.view{"BBMap: $it"}
  .concat(all_reports).flatten().collect()

  // Salmon 
  all_reports = ch_salmon_result.map{it -> it[1]}.ifEmpty([])
  //.view{"Salmon quant: $it"}
  .concat(all_reports).flatten().collect()

  // Kallisto
  all_reports = ch_kallisto_result.map{it -> it[2]}.ifEmpty([])
  //.view{"Kallisto quant: $it"}
  .concat(all_reports).flatten().collect()

  // Picard
  all_reports = ch_picard_rnametrics.map{it -> it[2]}.ifEmpty([])
  //.view{"Picard RNAseqmetrics: $it"}
  .concat(all_reports).flatten().collect()
  all_reports = ch_alignment_markdups.map{it -> it[4]}.ifEmpty([])
  //.view{"Picard markduplicates: $it"}
  .concat(all_reports).flatten().collect()
  
  // RNASEQC
  all_reports = ch_rnaseqc.map{it -> it[2]}.ifEmpty([])
  //.view{"RNA-SeQC: $it"}
  .concat(all_reports).flatten().collect()

  //bowtie2_err = ch_alignment_output.map{it -> it[2]}.collect().ifEmpty([])
  //kraken_err = ch_kraken2_output.map{it -> it[2]}.collect().ifEmpty([])
  //bracken_err = ch_bracken_output.map{it -> it[4]}.collect().ifEmpty([])
  //metaphlan_err = ch_metaphlan_output.map{it -> it[1]}.collect().ifEmpty([])
  //megahit_err = ch_megahit_output.map{it -> it[4]}.collect().ifEmpty([])
  //ch_metaquast_tsv = ch_metaquast_output.map{it -> it[3]}.collect().ifEmpty([])

  //Call multiQC process
  multiQC(params.multiQC.configyaml,
            all_reports
            )
  ch_multiqc_out = multiQC.out //.map{it -> it[1]}
emit:
ch_multiqc_out
}
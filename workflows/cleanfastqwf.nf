


include { doTrimmomatic } from '../modules/trimmomatic'
include { doCutadapt } from '../modules/cutadapt'
include { doBBduk } from '../modules/bbduk'
include { getFastQCIllumina } from '../modules/fastqc'


workflow CLEANFASTQ {
  take: ch_rawfastq
  main:
 //Trim reads
  ch_fastq_processed = ch_rawfastq
  if(params.workflows.do_trim && params.workflows.trimming_tool == "trimmomatic" ){ //With trimmomatic
    doTrimmomatic(ch_rawfastq)
    ch_fastq_processed = doTrimmomatic.out
     .view{ "Illumina trimmed reads with Trimmomatic: $it" }
     ch_fastq_processed_paired = ch_fastq_processed.map{it -> tuple(it[0], it[1])}
     .view{ "trimmed fastq paired only: $it" }
  }else if (params.workflows.do_trim && params.workflows.trimming_tool == "cutadapt" ){ //With cutadapt
     doCutadapt(ch_rawfastq)
    ch_fastq_processed = doCutadapt.out
     .view{ "Illumina trimmed reads with Cutadapt: $it" }
     ch_fastq_processed_paired = ch_fastq_processed.map{it -> tuple(it[0], it[1])}
     .view{ "trimmed fastq paired only: $it" } 
   
  }else if (params.workflows.do_trim && params.workflows.trimming_tool == "bbduk" ){ //With bbduk
     doBBduk(ch_rawfastq)
    ch_fastq_processed = doBBduk.out
     .view{ "Illumina trimmed reads with BBduk: $it" }
     ch_fastq_processed_paired = ch_fastq_processed.map{it -> tuple(it[0], it[1])}
     .view{ "trimmed fastq paired only: $it" } 
   
  }else{ //Do not trim
    ch_fastq_processed = ch_rawfastq
    ch_fastq_processed_paired = ch_rawfastq
    .view{ "input fastq into processed channel: $it" }
  }

  //Initialize channels 
  ch_fastqc = Channel.from([])
  ch_flatfastq = Channel.from([])

  //Get a flattened list of raw fastq
  if(params.getFastQCIllumina.do_fastqc_raw){
    ch_flatfastq = ch_rawfastq.map{it -> it[1]}.flatten()
  }

  //Add trimmed fastq to the flat fastq channel to perform FastQC
    if(params.getFastQCIllumina.do_fastqc_trim &&  params.workflows.do_trim){
      ch_flatfastq = ch_fastq_processed.map{it -> it[1]}
      .flatten().concat(ch_flatfastq)
    }
    if(params.getFastQCIllumina.do_fastqc_trim_single && params.workflows.do_trim && params.workflows.trimming_tool == "trimmomatic"){
      ch_flatfastq = ch_fastq_processed.map{it -> it[2]}
      .flatten().concat(ch_flatfastq)
    }

  // Do FastQC for all fastq files
  if(params.getFastQCIllumina.do_fastqc_raw | params.getFastQCIllumina.do_fastqc_trim){
    getFastQCIllumina(ch_flatfastq)
    //.view{ "getFastQCIllumina - FastQC reports: $it" }
    ch_fastqc = getFastQCIllumina.out
  }


  emit:
  ch_fastq_processed
  ch_fastq_processed_paired
  ch_fastqc

}
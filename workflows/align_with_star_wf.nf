include { alignSTAR } from '../modules/star-align'
include { buildindexSTAR } from '../modules/star-buildindex'
include { alignSTAR2ndPass } from '../modules/star-align-2ndpass'

workflow ALIGN_WITH_STAR {
    take: ch_fastq_processed_paired
    main:

    //Create STAR index
    if(params.buildindexSTAR.do_index){  
      buildindexSTAR(params.buildindexSTAR.index_name,
                    params.buildindexSTAR.fasta,
                    params.buildindexSTAR.annot)

      ch_star_index = buildindexSTAR.out
       //.view{ "STAR index created: $it" }
    }else{
      ch_star_index = params.alignSTAR.index
    }

    alignSTAR(ch_star_index, ch_fastq_processed_paired)

    ch_star_result = alignSTAR.out
        //.view{ "STAR full result: $it" }
    ch_star_bam = ch_star_result.map{it -> tuple("STAR", it[0], it[1], it[2])}
        //.view{ "STAR BAM only: $it" }
    ch_star_bam_bytranscript = ch_star_result.map{it -> tuple("STARt", it[0], it[3], it[4])}
        //.view{ "STAR BAM only by transcript: $it" }

    // Align with STAR 2nd pass, with annotated splice junctions from the 1st pass
    if(params.alignSTAR2ndPass.do_star){
        ch_star2ndpass_input_SJ = ch_star_result
            .map{it -> tuple(it[5])}
            .collect()
            //.view{ "STAR splice junctions collect: $it" }

        alignSTAR2ndPass(ch_star_index, ch_fastq_processed_paired, ch_star2ndpass_input_SJ)

        ch_star_2ndpass_result = alignSTAR2ndPass.out
            //.view{ "STAR 2nd PASS full result: $it" }
        ch_star_2ndpass_bam = ch_star_2ndpass_result.map{it -> tuple("STAR2", it[0], it[1], it[2])}
            .view{ "STAR 2nd PASS BAM only: $it" }
        ch_star_2ndpass_bam_bytranscript = ch_star_2ndpass_result.map{it -> tuple("STARt2", it[0], it[3], it[4])}
            .view{ "STAR 2nd PASS BAM only by transcript: $it" }
    } else{
        ch_star_2ndpass_result = Channel.from([])
        ch_star_2ndpass_bam = Channel.from([])
        ch_star_2ndpass_bam_bytranscript = Channel.from([])
    }// end STAR 2nd pass

  emit:
  ch_star_result
  ch_star_bam
  ch_star_bam_bytranscript
  ch_star_2ndpass_result
  ch_star_2ndpass_bam
  ch_star_2ndpass_bam_bytranscript
}
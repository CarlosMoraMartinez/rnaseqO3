include { buildindexHISAT2 } from '../modules/hisat2-buildindex'
include { alignHISAT2 } from '../modules/hisat2-align'
include { alignSTAR } from '../modules/star-align'
include { buildindexSTAR } from '../modules/star-buildindex'
include { alignSTAR2ndPass } from '../modules/star-align-2ndpass'
include { quantStringtie2 } from '../modules/stringtie2-quant'
include { mergeStringtie2 } from '../modules/stringtie2-merge'
include { quantSalmon } from '../modules/salmon-quant'
include { mergeSalmon } from '../modules/salmon-merge'
include { mapdecoySalmontools } from '../modules/salmontools-mapdecoy'
include { buildindexSalmon } from '../modules/salmon-buildindex'

workflow ALIGN {
  take: ch_fastq_processed_paired
  main:
  //Align using HISAT2
  ch_hisat2_result = Channel.from([])
  ch_hisat2_bam = Channel.from([])

  if(params.alignHISAT2.do_hisat2){ 
    if(params.buildindexHISAT2.do_index){
      buildindexHISAT2(
        params.buildindexHISAT2.index_dir,
        params.buildindexHISAT2.index_name,
        params.buildindexHISAT2.fasta,
        params.buildindexHISAT2.annot)
      ch_hisat2_index = buildindexHISAT2.out
        .map{it -> it[0]}
        //.view{ "HISAT2 index created: $it" }
    }else{
      ch_hisat2_index = params.alignHISAT2.index
    }
    alignHISAT2(ch_hisat2_index, ch_fastq_processed_paired)
    ch_hisat2_result = alignHISAT2.out
     //.view{ "HISAT2 full result: $it" }
    ch_hisat2_bam = ch_hisat2_result.map{it -> tuple("HISAT2", it[0], it[1], it[2])}
      //.view{ "HISAT2 BAM only: $it" }
  } // end HISAT2

   //Align using STAR
  ch_star_result = Channel.from([])
  ch_star_bam = Channel.from([])
  ch_star_2ndpass_result = Channel.from([])
  ch_star_2ndpass_bam = Channel.from([])

  if(params.alignSTAR.do_star){ 
    if(params.buildindexSTAR.do_index){
      
      buildindexSTAR(params.buildindexSTAR.index_name,
                    params.buildindexSTAR.fasta,
                    params.buildindexSTAR.annot)

      ch_star_index = buildindexSTAR.out
       //.view{ "STAR index created: $it" }
    }else{
      ch_star_index = params.alignSTAR.index
    }

    alignSTAR(ch_star_index, 
              ch_fastq_processed_paired)

     ch_star_result = alignSTAR.out
     //.view{ "STAR full result: $it" }
     ch_star_bam = ch_star_result.map{it -> tuple("STAR", it[0], it[1], it[2])}
     //.view{ "STAR BAM only: $it" }

    if(params.alignSTAR2ndPass.do_star){
    ch_star2ndpass_input_SJ = ch_star_result
        .map{it -> tuple(it[3])}
        .collect()
       //.view{ "STAR splice junctions collect: $it" }

    alignSTAR2ndPass(ch_star_index, 
                    ch_fastq_processed_paired, 
                    ch_star2ndpass_input_SJ)

    ch_star_2ndpass_result = alignSTAR2ndPass.out
    //.view{ "STAR 2nd PASS full result: $it" }
     ch_star_2ndpass_bam = ch_star_2ndpass_result.map{it -> tuple("STAR2", it[0], it[1], it[2])}
    //.view{ "STAR 2nd PASS BAM only: $it" }
    }// end STAR 2nd pass
  }
  // Merge alignments
  ch_alignment_all = ch_hisat2_bam
    .concat(ch_star_bam)
    .concat(ch_star_2ndpass_bam)
    //.view{ "All alignments concat: $it" }

  ch_stringtie_results = Channel.from([])
  if(params.quantStringtie2.do_stringtie){
    quantStringtie2(ch_alignment_all)
    ch_stringtie_results = quantStringtie2.out
    //.view{ "Stringtie results: $it" }
    ch_stringtie_results_grouped = ch_stringtie_results
    .map{it -> tuple(it[0], tuple(it[2]))}
    .groupTuple(by: 0)
    .map{it -> tuple(it[0], it[1].flatten())}
    //.view{ "Stringtie results grouped: $it" }

    mergeStringtie2(ch_stringtie_results_grouped)
    ch_stringtie_results_merged = mergeStringtie2.out
      //.view{ "Stringtie merge prepDE.py results: $it" }
  }
  
  if(params.quantSalmon.do_salmon){
    if(params.buildindexSalmon.do_index){
      if(params.buildindexSalmon.mode == "partial_decoy"){
        mapdecoySalmontools(
          params.mapdecoySalmontools.name,
          params.mapdecoySalmontools.genome_fasta,
          params.mapdecoySalmontools.transcripts_fasta,
          params.mapdecoySalmontools.annot
        )
        ch_salmon_partialdecoy = mapdecoySalmontools.out
          .view{ "Salmon partial decoy created: $it" }
        ch_partial_decoy_fasta = ch_salmon_partialdecoy.map{it -> it[0]}
          .view{ "Salmon partial decoy FASTA: $it" }
        ch_partial_decoy_txt = ch_salmon_partialdecoy.map{it -> it[1]}
          .view{ "Salmon partial decoy txt: $it" }

        buildindexSalmon(
          params.buildindexSalmon.name,
          params.buildindexSalmon.mode,
          params.buildindexSalmon.genome_fasta,
          ch_partial_decoy_fasta,
          ch_partial_decoy_txt
        )
        ch_salmon_index = buildindexSalmon.out
         .view{ "Salmon index created with partial decoy: $it" }

      }else{
        buildindexSalmon(
          params.buildindexSalmon.name,
          params.buildindexSalmon.mode,
          params.buildindexSalmon.genome_fasta,
          params.buildindexSalmon.transcripts_fasta,
          params.buildindexSalmon.annot //Not used for anything, decoy when mode=partial_decoy
        )
        ch_salmon_index = buildindexSalmon.out
         .view{ "Salmon index created: $it" }
      }//if partial decoy
    }else{
      ch_salmon_index = params.quantSalmon.index
    }// if build index
    quantSalmon(ch_salmon_index, ch_fastq_processed_paired)
    ch_salmon_result = quantSalmon.out
      .view{ "Salmon results: $it" }
    ch_salmon_names = ch_salmon_result.map{it -> it[0]}.collect().map{ it -> it.join(' ')}
    ch_salmon_dirs = ch_salmon_result.map{it -> it[1]}.collect()
    mergeSalmon("salmon_" + params.buildindexSalmon.mode,
                ch_salmon_names, 
                ch_salmon_dirs)
    ch_salmon_merged = mergeSalmon.out
      .view{ "Salmon results merged: $it" }

  } //if Do salmon

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
  ch_salmon_merged
}
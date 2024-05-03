include { ALIGN_WITH_HISAT2 } from './align_with_hisat2_wf'
include { alignHISAT2 } from '../modules/hisat2-align'
include { alignSTAR } from '../modules/star-align'
include { buildindexSTAR } from '../modules/star-buildindex'
include { alignSTAR2ndPass } from '../modules/star-align-2ndpass'
include { quantStringtie2 } from '../modules/stringtie2-quant'
include { mergeStringtie2 } from '../modules/stringtie2-merge'
include { quantSalmon } from '../modules/salmon-quant'
include { quantBamSalmon } from '../modules/salmon-quant-frombam'
include { mergeSalmon } from '../modules/salmon-merge'
include { mapdecoySalmontools } from '../modules/salmontools-mapdecoy'
include { buildindexSalmon } from '../modules/salmon-buildindex'

workflow ALIGN {
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
  ch_star_2ndpass_result = Channel.from([])
  ch_star_2ndpass_bam = Channel.from([])
  ch_star_bam_bytranscript = Channel.from([])
  ch_star_2ndpass_bam_bytranscript = Channel.from([])
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
     ch_star_bam_bytranscript = ch_star_result.map{it -> tuple("STARt", it[0], it[3], it[4])}
     //.view{ "STAR BAM only by transcript: $it" }

    if(params.alignSTAR2ndPass.do_star){
    ch_star2ndpass_input_SJ = ch_star_result
        .map{it -> tuple(it[5])}
        .collect()
       //.view{ "STAR splice junctions collect: $it" }

    alignSTAR2ndPass(ch_star_index, 
                    ch_fastq_processed_paired, 
                    ch_star2ndpass_input_SJ)

    ch_star_2ndpass_result = alignSTAR2ndPass.out
    //.view{ "STAR 2nd PASS full result: $it" }
     ch_star_2ndpass_bam = ch_star_2ndpass_result.map{it -> tuple("STAR2", it[0], it[1], it[2])}
    .view{ "STAR 2nd PASS BAM only: $it" }
    ch_star_2ndpass_bam_bytranscript = ch_star_2ndpass_result.map{it -> tuple("STARt2", it[0], it[3], it[4])}
    .view{ "STAR 2nd PASS BAM only by transcript: $it" }
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
  ch_salmon_result = Channel.from([])
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
          //.view{ "Salmon partial decoy created: $it" }
        ch_partial_decoy_fasta = ch_salmon_partialdecoy.map{it -> it[0]}
          //.view{ "Salmon partial decoy FASTA: $it" }
        ch_partial_decoy_txt = ch_salmon_partialdecoy.map{it -> it[1]}
          //.view{ "Salmon partial decoy txt: $it" }

        buildindexSalmon(
          params.buildindexSalmon.name,
          params.buildindexSalmon.mode,
          params.buildindexSalmon.genome_fasta,
          ch_partial_decoy_fasta,
          ch_partial_decoy_txt
        )
        ch_salmon_index = buildindexSalmon.out
         //.view{ "Salmon index created with partial decoy: $it" }

      }else{
        buildindexSalmon(
          params.buildindexSalmon.name,
          params.buildindexSalmon.mode,
          params.buildindexSalmon.genome_fasta,
          params.buildindexSalmon.transcripts_fasta,
          params.buildindexSalmon.annot //Not used for anything, decoy when mode=partial_decoy
        )
        ch_salmon_index = buildindexSalmon.out
         //.view{ "Salmon index created: $it" }
      }//if partial decoy
    }else{
      ch_salmon_index = params.quantSalmon.index
    }// if build index
    quantSalmon(ch_salmon_index, ch_fastq_processed_paired)
    ch_salmon_result = quantSalmon.out
      //.view{ "Salmon results: $it" }
  } //if Do salmon

  //Salmon quantify from bam
  ch_alignment_all_bytrans = Channel.from([])
  ch_salmon_aln_result = Channel.from([])
  ch_salmon_aln_names = Channel.from([])
  ch_salmon_aln_dirs = Channel.from([])
  if(params.alignSTAR.do_star && params.quantBamSalmon.do_salmon){
      ch_alignment_all_bytrans = ch_star_bam_bytranscript  //Concat STAR results from 1st and 2nd pass
            .concat(ch_star_2ndpass_bam_bytranscript)
      quantBamSalmon(ch_alignment_all_bytrans)  //Quantify with Salmon in alignment mode
      ch_salmon_aln_result = quantBamSalmon.out
      //.view{ "Salmon from Alignment results: $it" }
  }// Salmon quant from SAM file 

  //Merge all salmon results (salmon mapping mode + salmon alignment mode from STAR 1st and 2nd passes)
  ch_salmon_merged = Channel.from([])
  if((params.quantBamSalmon.do_salmon && params.alignSTAR.do_star) || params.quantSalmon.do_salmon){
      ch_salmon_aln_results_grouped = ch_salmon_result.map{it -> tuple("salmon-"+ params.buildindexSalmon.mode, it[0], it[1])}
        .concat(ch_salmon_aln_result)
        //.view{ "Salmon ALN concat: $it" }
        .multiMap{it ->
                      samplenames: tuple(it[0], tuple(it[1]))
                      salmonresults: tuple(it[0], tuple(it[2]))                 
                      }
        ch_salmon_aln_results_grouped.samplenames.groupTuple(by: 0).map{it2 -> tuple(it2[0], it2[1].flatten().join(' '))}
            .view{ "Salmon ALN names SEP: $it" }
            .set{salmon_grouped_samplenames}

        ch_salmon_aln_results_grouped.salmonresults.groupTuple(by: 0).map{it3 -> tuple(it3[0], it3[1].flatten())}
            .view{ "Salmon ALN results SEP: $it" }
            .set{salmon_grouped_resultdirs}
            
        ch_salmon_aln_results_grouped2 = salmon_grouped_samplenames
            .combine(salmon_grouped_resultdirs, by: 0)
            .view{ "Salmon ALN results MERGED: $it" }

        mergeSalmon(ch_salmon_aln_results_grouped2)
        ch_salmon_merged = mergeSalmon.out
        .view{ "Salmon results merged: $it" }
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
}
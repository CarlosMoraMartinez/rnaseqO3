include { quantSalmon } from '../modules/salmon-quant'
include { quantBamSalmon } from '../modules/salmon-quant-frombam'
include { mergeSalmon } from '../modules/salmon-merge'
include { mapdecoySalmontools } from '../modules/salmontools-mapdecoy'
include { buildindexSalmon } from '../modules/salmon-buildindex'

workflow QUANTIFY_WITH_SALMON {
  take: 
  ch_fastq_processed_paired
  ch_star_bam_bytranscript
  ch_star_2ndpass_bam_bytranscript
  main:

  // 1. Build Salmon index, if necessary. Otherwise use index passed as parameter 
  if(params.buildindexSalmon.do_index){
    // If the type of index is 'partial_decoy', use salmontools to calculate the partial decoy
    // Otherwise, build the index directly from transcriptome or transcriptome+genome fasta. 
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
  
  // 2. Quantify with Salmon directly from fastq, using index.
  ch_salmon_result = Channel.from([])
  if(params.quantSalmon.do_salmon){ 
    quantSalmon(ch_salmon_index, ch_fastq_processed_paired)
    ch_salmon_result = quantSalmon.out
    //.view{ "Salmon results: $it" }
  } //if Do salmon from fastq

  // 3. Quantify with Salmon from bam. Index not required.
  ch_salmon_aln_result = Channel.from([])
  if(params.workflows.do_star && params.quantBamSalmon.do_salmon){
    ch_alignment_all_bytrans = ch_star_bam_bytranscript  //Concat STAR results from 1st and 2nd pass
      .concat(ch_star_2ndpass_bam_bytranscript)
    quantBamSalmon(ch_alignment_all_bytrans)  //Quantify with Salmon in alignment mode
    ch_salmon_aln_result = quantBamSalmon.out
    //.view{ "Salmon from Alignment results: $it" }
  }// Salmon quant from SAM file 

  // 4. Merge all salmon results (salmon mapping mode + salmon alignment mode from STAR 1st and 2nd passes)
  ch_salmon_merged = Channel.from([])
  if((params.quantBamSalmon.do_salmon && params.alignSTAR.do_star) || params.quantSalmon.do_salmon){
    // First concatenate results from salmon mapping and alignment mode,
    // and use multimap to separate into a channel with sample names and another with results
    ch_salmon_aln_results_grouped = ch_salmon_result
      .map{it -> tuple("salmon-"+ params.buildindexSalmon.mode, it[0], it[1])}
      .concat(ch_salmon_aln_result)
      .multiMap{it ->
            samplenames: tuple(it[0], tuple(it[1]))
            salmonresults: tuple(it[0], tuple(it[2]))         
            }
    // Group the names channel by process of origin, and join names into a string separated by spaces
    // result is a channel of: [process_type, "name1 name2 ..."]
    ch_salmon_aln_results_grouped.samplenames.groupTuple(by: 0)
      .map{it2 -> tuple(it2[0], it2[1].flatten().join(' '))}
      //.view{ "Salmon ALN names SEP: $it" }
      .set{salmon_grouped_samplenames}

    // Group the results channel by process of origin into a list
    // result is a channel of: [process_type, [dir1, dir2, ...]]
    ch_salmon_aln_results_grouped.salmonresults.groupTuple(by: 0)
      .map{it3 -> tuple(it3[0], it3[1].flatten())}
      //.view{ "Salmon ALN results SEP: $it" }
      .set{salmon_grouped_resultdirs}
    
    // Combine the previous two channels by process type
    // result is a channel of: [process_type, "name1 name2 ...", [dir1, dir2, ...]]
    ch_salmon_aln_results_grouped2 = salmon_grouped_samplenames
      .combine(salmon_grouped_resultdirs, by: 0)
      //.view{ "Salmon ALN results COMBINED: $it" }

    // Finally call the process to merge salmon results of different samples, by process type
    mergeSalmon(ch_salmon_aln_results_grouped2)
    ch_salmon_merged = mergeSalmon.out
    //.view{ "Salmon merge results: $it" }
  }

  emit:
  ch_salmon_result
  ch_salmon_aln_result
  ch_salmon_merged
}
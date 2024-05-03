include { buildindexKallisto } from '../modules/kallisto-buildindex'
include { quantKallisto } from '../modules/kallisto-quant'

workflow QUANTIFY_WITH_KALLISTO {
  take: ch_fastq_processed_paired

  main:

  if(params.buildindexKallisto.do_index){
    buildindexKallisto(
      params.buildindexKallisto.index_name,
      params.buildindexKallisto.transcripts_fasta
      )
    ch_kallisto_index = buildindexKallisto.out
      //.view{ "KALLISTO index created: $it" }
  }else{
    ch_kallisto_index = params.quantKallisto.index
  }
  quantKallisto(ch_kallisto_index, ch_fastq_processed_paired)
  ch_kallisto_result = quantKallisto.out
    //.view{ "KALLISTO full result: $it" }

  emit:
  ch_kallisto_result
}
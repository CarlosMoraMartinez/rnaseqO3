
process buildRefflatFromGTF{
  label 'mg00_build_refflat'
  conda params.buildRefflatFromGTF.conda
  cpus params.resources.buildRefflatFromGTF.cpus
  memory params.resources.buildRefflatFromGTF.mem
  queue params.resources.buildRefflatFromGTF.queue
  clusterOptions params.resources.buildRefflatFromGTF.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg00_build_refflat", mode: 'symlink'
  input:
  path(genome_fasta)
  path(annot)
  
  output:
  path('*.refflat')

  shell:
  '''
  outname=$(basename -s .gtf !{annot})'.refflat'
        
  # It fails when symbolic links are used, therefore, use val inputs with gtftorefflat
  # Actually, doesnt seem to work most of the times. 
  #java -jar !{params.software.gtftorefflat_path} \
  #-g !{annot} \
  #-R !{genome_fasta} \
  #-r $outname 2>gtftorefflat.err
  
  gtfToGenePred !{annot} tmp.txt 2>gtftorefflat.err
  
  awk '
    BEGIN { FS=OFS="\t" } 
    {
        $11 = $1; 
        print $1, $11, $2, $3, $4, $5, $6, $7, $8, $9, $10 
    }' tmp.txt | 
    awk '
    BEGIN { FS=OFS="\t" } 
    !($3 ~ /chrNT/)
    ' > $outname
  rm tmp.txt
  '''
}


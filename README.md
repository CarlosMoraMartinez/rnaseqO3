# RNA-seqO3
## An RNA-seq pipeline to quantify expression from short reads

This pipeline combines some of the most commond RNA-seq analysis tools. Currently, only one trimming tool can be used, but all aligners and quantifiers are combined (i.e., if using 2 aligners and 2 quantifiers, you will end up with 4 count matrices). 

Currently, the following trimming tools are available:

- Trimmomatic
- BBDuk
- Cutadapt

The following software combinations are available:

- HISAT2 + Stringtie2
- HISAT2 + featureCounts
- HISAT2 + HTSeq
- STAR 1 pass + Stringtie2
- STAR 1 pass + featureCounts
- STAR 1 pass + HTSeq
- STAR 1 pass + Salmon in alignment mode
- STAR 2 passes + Stringtie2
- STAR 2 passes + featureCounts
- STAR 2 passes + HTSeq
- subread + Stringtie2
- subread + featureCounts
- subread + HTSeq
- subread + Salmon in alignment mode
- Kallisto
- Salmon

All tools can be configured by changing the config file. There is support to build the three index types in salmon (from transcriptome without decoy, with full decoy and with partial decoy, by calling the [generateDecoyTranscriptome.sh](https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh) script internally). HTSeq can be run with combinations of different options in the same run.  

# Running the pipeline

The pipeline is prepared to run locally or in a Slurm cluster. You can switch between configurations by using different config files (in the **/config** directory)

```
nextflow run main.nf -c config/run_samples_cluster.config -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.html
```

# Dependencies

The pipeline uses different Docker containers for each step; some are pre-build and some are built here. All the dockerfiles can be found in the **docker/** directory. The commands to build them are in **docker/build_image.sh**, and the commands to download the rest are in **docker/docker_images.sh**. By changing the configuration file it is easy to use conda instead of Docker.

# Structure of the repository

Individual processes (e.g., call FastQC, call HISAT2, etc) are in individual files in the **modules/** directory.
Workflows use a set of related processes (e.g., index transcriptome -> call Salmon). Each workflow is in an individual file in the **workflows/** directory. Some workflows call other workflows.

Configuration files are in **config/**, and **sbatch** files to launch the pipeline in a server are in **sbatch/**. The **scripts/** directory contains some useful extra scripts. 
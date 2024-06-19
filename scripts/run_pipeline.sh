sudo nextflow run main.nf -c config/run_samples_local_cmm_docker.config -resume

sudo nextflow run main.nf -c config/run_samples_local_cmm_docker.config -profile conda -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.html
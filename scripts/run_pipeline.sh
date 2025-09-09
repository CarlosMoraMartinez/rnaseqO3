sudo nextflow run main.nf -c config/run_samples_local_cmm_docker.config -resume

sudo nextflow run main.nf -c config/run_samples_local_cmm_docker.config -profile conda -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.html

sudo /home/servidor/PROGRAMAS/nextflow run main.nf -c config/run_mm1_local_1.config -profile conda -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.html


sudo rm *.html
sudo /home/servidor/PROGRAMAS/nextflow run main.nf -c config/run_mm1_local_sim.config -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.htmlne_dag.html

# git switch -c feature/evaMouseSept2025_2 origin/feature/evaMouseSept2025_2
# /home/servidor/bin/nextflow on server 3

sudo /home/servidor/bin/nextflow  run main.nf -c config/run_evaMouseSept25_docker.config -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.htmlne_dag.html
# dedup + STAR
nextflow run main.nf --dedup -w scratch -params-file params_pe.yaml -profile singularity
mv results ../dedup_star_results
rm -r scratch
# dedup + bowtie2
nextflow run main.nf --bowtie2 --dedup -w scratch -params-file params_bowtie2.yaml -profile singularity
mv results ../dedup_bowtie_results
rm -r scratch
# bowtie 2
nextflow run main.nf --bowtie2 -w scratch -params-file params_bowtie2.yaml -profile singularity
mv results ../bowtie_results
rm -r scratch
# STAR 
nextflow run main.nf -w scratch -params-file params_pe.yaml -profile singularity
mv results ../star_results
rm -r scratch

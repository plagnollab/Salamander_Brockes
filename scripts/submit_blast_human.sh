#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -l tmem=1.9G,h_vmem=1.9G
#$ -pe smp 1
#$ -l h_rt=20:0:0
#$ -V
#$ -R y
#$ -t 1-158
#$ -tc 80


fasta=/SAN/biomed/biomed5/biomed5/salamander_paper/split_genome/AxoJuveBlastema_RSEM-EVAL-1.4_Trinity_${SGE_TASK_ID}.fasta
blastx=/cluster/project8/vyp/vincent/Software/ncbi-blast-2.2.29+/bin/blastx
database=/scratch2/vyp-scratch2/reference_datasets/blast/human_protein/human_prot

if [ ! -s data/blast/blast_output_${SGE_TASK_ID}.tab ]; then
    $blastx -db $database -query $fasta -out data/blast/blast_output_${SGE_TASK_ID}.tab -outfmt 6 -num_threads 1
fi

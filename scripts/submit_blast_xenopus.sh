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

##SGE_TASK_ID=1

fasta=/SAN/biomed/biomed5/biomed5/salamander_paper/split_genome/AxoJuveBlastema_RSEM-EVAL-1.4_Trinity_${SGE_TASK_ID}.fasta
blastx=/cluster/project8/vyp/vincent/Software/ncbi-blast-2.2.29+/bin/blastx
database=protein_databases/xenopus/xenopus

if [ ! -s data/blast_xenopus/blast_output_${SGE_TASK_ID}.tab ]; then
    $blastx -db $database -query $fasta -out data/blast_xenopus/blast_output_${SGE_TASK_ID}.tab -outfmt 6 -num_threads 1
fi



# Salamander Brockes project

## Experimental design:

* 15 samples
* 3 time points (D0, D15, D150)


Species: axolotl

Reference species used for annotation of genes:
* zebrafish
* xenopus
* human


## Overview of analysis

Read aligned to reference contigs using RSEM.

Overview of raw counts can be found here:
~~~~bash
zless -S data/Genecounts_final_table.txt.gz 
~~~~




## List of scripts

Run the deseq2 analysis with the following script:
~~~~bash
R < scripts/deseq2_analysis.R --no-save
~~~~

Annotate the contig with the pythong script
~~~~bash
python scripts/annotate_data_v2_all_species.py
~~~~

Final clean up of the data using
~~~~bash
R < scripts/clean_up_diff_expression.R --no-save
~~~~
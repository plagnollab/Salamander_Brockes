for i in `seq 1 159`; do if [ ! -s  data/blast_human/blast_output_${i}.tab ]; then echo "Missing data/blast_human/blast_output_${i}.tab";fi; done


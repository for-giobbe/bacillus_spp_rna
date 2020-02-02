for i in $(ls -l /scratch/rh35/gf6030/snake/raw_reads/BAT*.gz | awk -F "/" '{print $NF}' | awk -F "_" '{print $1"_"$2"_"$3"_"$4}' | sort -u);

do species=$(echo $i | awk -F "_" '{print $1}')

/scratch/rh35/gf6030/Miniconda3_x86_64/opt/trinity-2.4.0/util/align_and_estimate_abundance.pl --seqType fq --left "/scratch/rh35/gf6030/snake/raw_reads/"$i"_1.fastq.gz" --right "/scratch/rh35/gf6030/snake/raw_reads/"$i"_2.fastq.gz" --transcripts "/scratch/rh35/gf6030/snake/expression_quantification/"$species".def.nuc.fa" --est_method RSEM  --aln_method bowtie2 --trinity_mode  --output_dir "/scratch/rh35/gf6030/snake/expression_quantification/"$i --prep_reference

done

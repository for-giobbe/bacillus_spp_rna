
	for i in ../alignemnts/5_sp_alignements/*fasta; 

		do printf "$i "; 

		for j in $(grep ">" $i);

			do transcript=$(echo $j | awk -F "_i" '{print $1}' | tr -d ">"); 

			species=$(echo $j | awk -F ".p1|.p2" '{print $2}' | tr -d "_1"); 

			t=$(grep $transcript DE_DESeq2_$species/*results | awk '{print $1"\t"$6"\t"$NF}'); 

			printf "$species $t "; 

		done ; 

	echo -e "\n"; 

	done


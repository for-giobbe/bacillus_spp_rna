	for i in *aln;

		do OG=$(echo $i | awk -F "." '{print $1}'); printf $OG;

		for j in $(grep ">" $i);

			do transcript=$(echo $j | awk -F "_i" '{print $1}' | tr -d ">");

			species=$(echo $j | awk -F ".p[0-9]" '{print $2}' | tr -d "_1");

			t=$(grep $transcript $species".result" | awk '{print $1" "$7" "$NF}');

			printf " $species  $t" >> test.tmp;

		done ;

	printf " $OG \n" >> test.tmp; echo -e "";

	done;



	while read line; 

		do echo $line | awk '/BAT TR*/ && /BRO TR*/ && /BGM TR*/'; 

	done < test.tmp > orthogroups-de_wout_m.tab;




	while read line; do transcript=$(echo $line | awk '{print $10}'); missing_line=$(grep $transcript BGM_m_only.results | awk '{print $1" "$7" "$NF}'); echo BGM_M $missing_line $line >> orthogroups-de_with_m.tmp; done < orthogroups-de_wout_m.tab

	while read line; do echo $line | awk '/BGM_M TR*/ && /BAT TR*/ && /BRO TR*/ && /BGM TR*/'; done < orthogroups-de_with_m.tmp > orthogroups-de_with_m.tab


	echo "BAT_transcript BAT_logFC BAT_padj BRO_transcript BRO_logFC BRO_padj BGM_transcript BGM_logFC BGM_padj OG" > orthogroups-de_wout_m.ref.tab
	cat orthogroups-de_wout_m.tab | awk '{print $2" "$3" "$4" "$6" "$7" " $8" "$10" "$11" "$12" "$13}' >> orthogroups-de_wout_m.ref.tab
        echo "BAT_transcript BAT_logFC BAT_padj BRO_transcript BRO_logFC BRO_padj BGM_F_transcript BGM_F_logFC BGM_F_padj  BGM_M_transcript BGM_M_logFC BGM_M_padj OG" > orthogroups-de_with_m.ref.tab
	cat orthogroups-de_with_m.tab | awk '{print $6" "$7" " $8" "$10" "$11" "$12" "$14" "$15" "$16" "$2" "$3" "$4" "$17}' >> orthogroups-de_with_m.ref.tab

        rm *tmp

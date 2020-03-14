while read line;

	do

	OG=$(echo $line | awk '{print $NF}');

	BAT_transcript=$(echo $line | awk '{print $1}');
	BRO_transcript=$(echo $line | awk '{print $4}');
	BGM_transcript=$(echo $line | awk '{print $7}');

	BGM_tmp_GO=$(grep $BGM_transcript 3sp_orthologs/BGM_argot_3sp.tsv | awk '{print $2}' | awk 'BEGIN { ORS = " " } { print }' | sed 's/ /, /g') ;
        BRO_tmp_GO=$(grep $BRO_transcript 3sp_orthologs/BRO_argot_3sp.tsv | awk '{print $2}' | awk 'BEGIN { ORS = " " } { print }' | sed 's/ /, /g') ;
        BAT_tmp_GO=$(grep $BAT_transcript 3sp_orthologs/BAT_argot_3sp.tsv | awk '{print $2}' | awk 'BEGIN { ORS = " " } { print }' | sed 's/ /, /g') ;

	echo -e "\n" >> 3sp_OG2GO.tab
	echo -ne "$OG\t $BGM_tmp_GO, $BRO_tmp_GO, $BAT_tmp_GO" | sed 's/, $//g' | sed 's/ , //g' >> 3sp_OG2GO.tab

	done < orthogroups-de_with_m.ref.tab

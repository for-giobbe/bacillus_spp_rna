echo -e "OG BAT_1+2_gc BAT_3rd_gc BRO_1+2_gc BRO_3rd_gc BGM_1+2_gc BGM_3rd_gc" > 4sp.gc3.reformat.tab
while read line;

	do

	og=$(echo $line | awk '{print $NF}');

	echo $og

	BGM_transcript=$(echo $line | awk '{print $2}'); BGM_gc=$(grep $BGM_transcript BGM.ref.gc | awk '{print $2" "$3}');
	echo BGM $BGM_transcript $BGM_gc;

	BAT_transcript=$(echo $line | awk '{print $6}'); BAT_gc=$(grep $BAT_transcript BAT.ref.gc | awk '{print $2" "$3}');
	echo BAT $BAT_transcript $BAT_gc;

	BRO_transcript=$(echo $line | awk '{print $10}');BRO_gc=$(grep $BRO_transcript BRO.ref.gc | awk '{print $2" "$3}');
	echo BRO $BRO_transcript $BRO_gc;

	echo -e "$og $BAT_gc $BRO_gc $BGM_gc" >> 4sp.gc3.reformat.tab

	done < ../4sp_ortho-de/orthogroups-de_with_m.tab


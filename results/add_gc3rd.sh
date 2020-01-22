while read line; 
 do 
 BGM_transcript=$(echo $line | awk '{print $2}'); BGMGC=$(grep $BGM_transcript gc3/BGM.ref.gc); echo $BGMGC;
 BAT_transcript=$(echo $line | awk '{print $6}'); BATGC=$(grep $BAT_transcript gc3/BAT.ref.gc); echo $BATGC;
 BRO_transcript=$(echo $line | awk '{print $10}');BROGC=$(grep $BRO_transcript gc3/BRO.ref.gc); echo $BROGC;
 echo -e "$line $BATGC $BROGC $BGMGC" >> 3sp_result.tmp
 done < 3sp_ortho-de/orthogroups-de_with_m.tab


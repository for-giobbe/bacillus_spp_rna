while read line;

do

og=$(echo $line | awk '{print $NF}');

BGM_transcript=$(echo $line | awk '{print $2}'); BGM_het=$(grep $BGM_transcript BAT_var.tab); 
if [ -z "$BGM_het" ]; then BGM_het=NA; BGM_transcript=NA; fi; echo BGM $BGM_het;

BRO_transcript=$(echo $line | awk '{print $2}'); BRO_het=$(grep $BRO_transcript BRO_var.tab); 
if [ -z "$BRO_het" ]; then BRO_het=NA; BRO_transcript=NA; fi; echo BRO $BRO_het;

BAT_transcript=$(echo $line | awk '{print $2}'); BAT_het=$(grep $BAT_transcript BAT_var.tab); 
if [ -z "$BAT_het" ]; then BAT_het=NA; BAT_transcript=NA; fi; echo BAT $BAT_het;

#echo -e "$og $BAT_thre $BGM_het $BRO_het $BAT_het " >> het.reformat.tab

done < ../3sp_ortho-de/orthogroups-de_with_m.tab

echo -e "OG BAT_3_utr_gc BAT_5_utr_gc BRO_3_utr_gc BRO_5_utr_gc BGM_3_utr_gc BGM_5_utr_gc" > 4sp.utr.reformat.tab
while read line;
do
og=$(echo $line | awk '{print $NF}');
BGM_transcript=$(echo $line | awk '{print $2}'); BGM_thre=$(grep $BGM_transcript BGM.thre.utr.gc | awk '{print $2}'); BGM_five=$(grep $BGM_transcript BGM.five.utr.gc | awk '{print $2}');
if [ -z "$BGM_thre" ]; then BGM_thre=NA; BGM_transcript=NA; fi ;if [ -z "$BGM_five" ]; then BGM_five=NA; fi; #echo BGM $BGM_thre $BGM_five;
BAT_transcript=$(echo $line | awk '{print $6}'); BAT_thre=$(grep $BAT_transcript BAT.thre.utr.gc | awk '{print $2}'); BAT_five=$(grep $BAT_transcript BAT.five.utr.gc | awk '{print $2}');
if [ -z "$BAT_thre" ]; then BAT_thre=NA; BAT_transcript=NA; fi ;if [ -z "$BAT_five" ]; then BAT_five=NA; fi; #echo BAT $BAT_thre $BAT_five;
BRO_transcript=$(echo $line | awk '{print $10}');BRO_thre=$(grep $BRO_transcript BRO.thre.utr.gc | awk '{print $2}'); BRO_five=$(grep $BRO_transcript BRO.five.utr.gc | awk '{print $2}');
if [ -z "$BRO_thre" ]; then BRO_thre=NA; BRO_transcript=NA; fi ;if [ -z "$BRO_five" ]; then BRO_five=NA; fi; #echo BRO $BRO_thre $BRO_five;
echo -e "$og $BAT_thre $BAT_five $BRO_thre $BRO_five $BGM_thre $BGM_five" >> 4sp.utr.reformat.tab
done < ../4sp_ortho-de/orthogroups-de_with_m.tab


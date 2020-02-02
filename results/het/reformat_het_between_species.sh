echo -e "OG BAT_len BAT_nvar BGM_len BGM_nvar BRO_len BRO_nvar" > het.BGM_F.tab

while read line;

do

og=$(echo $line | awk '{print $NF}'); echo $og

BGM_transcript=$(echo $line | awk '{print $2}'); BGM_het=$(grep $BGM_transcript BGM_F_var.tab | awk '{print $2" "$3}');
if [ -z "$BGM_het" ]; then BGM_het="- 0"; BGM_transcript=NA; fi; echo BGM $BGM_het;

BRO_transcript=$(echo $line | awk '{print $6}'); BRO_het=$(grep $BRO_transcript BRO_F_var.tab | awk '{print $2" "$3}');
if [ -z "$BRO_het" ]; then BRO_het="- 0"; BRO_transcript=NA; fi; echo BRO $BRO_het;

BAT_transcript=$(echo $line | awk '{print $10}'); BAT_het=$(grep $BAT_transcript BAT_F_var.tab | awk '{print $2" "$3}');
if [ -z "$BAT_het" ]; then BAT_het="- 0"; BAT_transcript=NA; fi; echo BAT $BAT_het;

echo -e "$og $BAT_het $BGM_het $BRO_het" >> het.BGM_F.tab

done < ../3sp_ortho-de/orthogroups-de_with_m.tab

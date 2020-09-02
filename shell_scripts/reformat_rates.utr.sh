while read line;

do

og=$(echo $line | awk '{print $1}');

BGM_transcript=$(echo $line | awk '{print $5}'); 
BGM_thre=$(grep -A 1 $BGM_transcript BGM.utr.3.fasta | tail -1); 
BGM_five=$(grep -A 1 $BGM_transcript BGM.utr.5.fasta | tail -1);
if [ -z "$BGM_thre" ]; then : ; else echo -e ">BGM_$BGM_transcript \n $BGM_thre" >> $og.utr.3.fasta; fi;
if [ -z "$BGM_five" ]; then : ; else echo -e ">BGM_$BGM_transcript \n $BGM_five" >> $og.utr.5.fasta; fi;

BAT_transcript=$(echo $line | awk '{print $4}'); 
BAT_thre=$(grep -A 1 $BAT_transcript BAT.utr.3.fasta | tail -1); 
BAT_five=$(grep -A 1 $BAT_transcript BAT.utr.5.fasta | tail -1);
if [ -z "$BAT_thre" ]; then : ; else echo -e ">BAT_$BAT_transcript \n $BAT_thre" >> $og.utr.3.fasta; fi;
if [ -z "$BAT_five" ]; then : ; else echo -e ">BAT_$BAT_transcript \n $BAT_five" >> $og.utr.5.fasta; fi; 

BRO_transcript=$(echo $line | awk '{print $3}');
BRO_thre=$(grep -A 1 $BRO_transcript BRO.utr.3.fasta | tail -1);
BRO_five=$(grep -A 1 $BRO_transcript BRO.utr.5.fasta | tail -1);
if [ -z "$BRO_thre" ]; then : ; else echo -e ">BRO_$BRO_transcript \n $BRO_thre" >> $og.utr.3.fasta; fi;
if [ -z "$BRO_five" ]; then : ; else echo -e ">BRO_$BRO_transcript \n $BRO_five" >> $og.utr.5.fasta; fi;

PHY_transcript=$(echo $line | awk '{print $2}');
PHY_thre=$(grep -A 1 $PHY_transcript PHY.utr.3.fasta | tail -1);
PHY_five=$(grep -A 1 $PHY_transcript PHY.utr.5.fasta | tail -1);
if [ -z "$PHY_thre" ]; then : ; else echo -e ">PHY_$PHY_transcript \n $PHY_thre" >> $og.utr.3.fasta; fi;
if [ -z "$PHY_five" ]; then : ; else echo -e ">PHY_$PHY_transcript \n $PHY_five" >> $og.utr.5.fasta; fi;

echo $BGM_transcript $BAT_transcript $BRO_transcript $PHY_transcript

mafft $og.utr.3.fasta > $og.utr.3.mafft.fasta
Gblocks $og.utr.3.mafft.fasta -t=d
cat $og.utr.3.mafft.fasta-gb | tr -d " " | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' > $og.utr.3.mafft.ref.fasta

mafft $og.utr.5.fasta > $og.utr.5.mafft.fasta
Gblocks $og.utr.5.mafft.fasta -t=d
cat $og.utr.5.mafft.fasta-gb | tr -d " " | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' > $og.utr.5.mafft.ref.fasta


done < OG.tab

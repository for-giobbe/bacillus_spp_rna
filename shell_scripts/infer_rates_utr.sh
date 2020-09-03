echo OG PHY BAT BRO BGM > rates.utr.tab

for i in *fasta;

do

OG=$(echo $i | awk -F "." '{print $1}')

raxmlHPC -m GTRGAMMA -s $i -n $i -f e -t sp.tre

PHY=$(awk -F "PHY" '{print $2}' "RAxML_result."$i | awk -F ":" '{print $2}' | awk -F "," '{print $1}' | tr -d "()")

BAT=$(awk -F "BAT" '{print $2}' "RAxML_result."$i | awk -F ":" '{print $2}' | awk -F "," '{print $1}' | tr -d "()")

BRO=$(awk -F "BRO" '{print $2}' "RAxML_result."$i | awk -F ":" '{print $2}' | awk -F "," '{print $1}' | tr -d "()")

BGM=$(awk -F "BGM" '{print $2}' "RAxML_result."$i | awk -F ":" '{print $2}' | awk -F "," '{print $1}' | tr -d "()")

echo $OG $PHY $BAT $BRO $BGM >> rates.utr.tab

done

rm RAxML*

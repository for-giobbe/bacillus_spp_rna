echo "OG PHY BRO BAT BGM" > OG.tab; 

for i in *aln; do 

OG=$(echo $i | awk -F "." '{print $1}'); 

PHY=$(grep PHY $i | tr -d ">" | awk -F "_" '{print $1}'); 
BRO=$(grep BRO $i | tr -d ">" | awk -F "_" 'NF{NF-=1};1' | sed 's/ /_/g'); 
BAT=$(grep BAT $i | tr -d ">" | awk -F "_" 'NF{NF-=1};1' | sed 's/ /_/g'); 
BGM=$(grep BGM $i | tr -d ">" | awk -F "_" 'NF{NF-=1};1' | sed 's/ /_/g'); 

echo $OG $PHY $BRO $BAT $BGM >> OG.tab; 

done


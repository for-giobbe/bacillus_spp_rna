#for i in *.mafft.n.aln; do distmat -sequence $i -nucmethod 1 -outfile $i".distances"; done

echo "OG BRO-BGM BAT-BGM BAT-BRO PHY-BAT PHY-BRO PHY-BGM" > 4sp.distances.tab

for i in *distances;

	do

	OG=$(echo $i | awk -F "." '{print $1}')

	BRO_BGM=$(sed '10q;d' $i | awk '{print $2}')
	BAT_BGM=$(sed '9q;d' $i | awk '{print $3}')
	BAT_BRO=$(sed '9q;d' $i | awk '{print $2}')

    	PHY_BAT=$(sed '9q;d' $i | awk '{print $4}')
        PHY_BRO=$(sed '10q;d' $i | awk '{print $3}')
        PHY_BGM=$(sed '11q;d' $i | awk '{print $2}')

	echo "$OG $BRO_BGM $BAT_BGM $BAT_BRO $PHY_BAT $PHY_BRO $PHY_BGM" >> 4sp.distances.tab

done;

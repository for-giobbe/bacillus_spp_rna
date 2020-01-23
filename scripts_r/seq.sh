echo "OG BAT_W-S BAT_S-W BRO_W-S BRO_S-W BGM_W-S BGM_S-W BGC_BAT BGC_BRO BGC_BGM ENC_BAT ENC_BRO ENC_BGM SCUO_BAT SCUO_BRO SCUO_BGM" > substiution_mapping.tab

for i in *aln;

	do
	OG=$( echo "$i " | awk -F "." '{print $1}')
	echo -n "$OG " >> substiution_mapping.tab
	sed "s/substitute_this/$i/g" seq.r > r.tmp;
	Rscript r.tmp;
	cat result.tmp >> substiution_mapping.tab
	done

rm *tmp

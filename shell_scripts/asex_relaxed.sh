# first argument of this script is the cutoff difference of dNdS to considered a selection regime to be relaxed

grep "alternative" 4sp_pur_sel.tab > var.tab.tmp;

	while read line;

	do

	BGM_dNdS=$(echo $line | awk '{print $3}');

	BAT_dNdS=$(echo $line | awk '{print $8}');

	BRO_dNdS=$(echo $line | awk '{print $13}');

	diff_BGM_BAT=$(echo "$BAT_dNdS - $BGM_dNdS" | bc);

        diff_BGM_BRO=$(echo "$BRO_dNdS - $BGM_dNdS" | bc);

	if (( $(echo "$diff_BGM_BAT > $1 " | bc -l) )); then echo $line | awk '{print $1}' >> BAT.relax.lst; fi;

        if (( $(echo "$diff_BGM_BRO > $1 " | bc -l) )); then echo $line | awk '{print $1}' >> BRO.relax.lst; fi;

	if (( $(echo "$diff_BGM_BRO > $1 " | bc -l) && $(echo "$diff_BGM_BAT > 0.1 " | bc -l) )); then echo $line | awk '{print $1}' >> automictics.relax.lst; fi;

	done < var.tab.tmp

	rm *tmp

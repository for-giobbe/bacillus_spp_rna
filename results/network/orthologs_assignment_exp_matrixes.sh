while read line;

	do

	transcript=$(echo $line | awk '{print $1}');
	OG=$(awk -v var="$transcript" '($10==var)' orthogroups-de_with_m.ref.tab | awk '{print $NF}');
		if [ ! -z "$OG" ];
		then sed -i "s/$transcript/$OG/" BGM_RSEM_mf.TMM.EXPR.matrix.ref;
		fi;
	echo $transcript $OG;

done < BGM_RSEM_mf.TMM.EXPR.matrix.ref

rm *lst

while read line;

do

OG=$(echo $line | awk '{print $NF}'); echo $OG

BGM_F_pval=$(echo $line | awk '{print $9}'); bBGM_F_pval=$( printf '%.40f' $BGM_F_pval);

BGM_M_pval=$(echo $line | awk '{print $12}'); bBGM_M_pval=$( printf '%.40f' $BGM_M_pval);

BRO_pval=$(echo $line | awk '{print $6}'); bBRO_pval=$( printf '%.40f' $BRO_pval);

BAT_pval=$(echo $line | awk '{print $3}'); bBAT_pval=$( printf '%.40f' $BAT_pval);

	if (( $(echo "$bBGM_F_pval > 0.005" | bc -l) )) && (( $(echo "$bBGM_M_pval < 0.005" | bc -l) ));
	then echo $OG >> BGM_M_gonad.specific.genes.lst;
	fi;

        if (( $(echo "$bBGM_M_pval > 0.005" | bc -l) )) && (( $(echo "$bBGM_F_pval < 0.005" | bc -l) ));
        then echo $OG >> BGM_F_gonad.specific.genes.lst;
        fi;

        if (( $(echo "$bBGM_F_pval > 0.005" | bc -l) )) && (( $(echo "$bBRO_pval < 0.005" | bc -l) )) && (( $(echo "$bBAT_pval < 0.005" | bc -l) ));
        then echo $OG >> DE.in.BRO+BAT_no.BGM.lst; echo $OG
        fi;

        if (( $(echo "$bBGM_F_pval > 0.005" | bc -l) )) && (( $(echo "$bBRO_pval < 0.005" | bc -l) )) && (( $(echo "$bBAT_pval > 0.005" | bc -l) ));
        then echo $OG >> DE.in.BRO_no.BGM+BAT.lst; echo $OG
        fi;

        if (( $(echo "$bBGM_F_pval > 0.005" | bc -l) )) && (( $(echo "$bBRO_pval > 0.005" | bc -l) )) && (( $(echo "$bBAT_pval < 0.005" | bc -l) ));
        then echo $OG >> DE.in.BAT_no.BGM+BRO.lst; echo $OG
        fi;

        if (( $(echo "$bBGM_F_pval < 0.005" | bc -l) )) && (( $(echo "$bBRO_pval > 0.005" | bc -l) )) && (( $(echo "$bBAT_pval > 0.005" | bc -l) ));
        then echo $OG >> DE.in.BGM_no.BAT+BRO.lst; echo $OG
        fi;

	if (( $(echo "$bBGM_F_pval < 0.005" | bc -l) )) && (( $(echo "$bBAT_pval > 0.005" | bc -l) ));
        then echo $OG >> DE.in.BGM_no.BAT.lst; echo $OG
        fi;

	if (( $(echo "$bBGM_F_pval < 0.005" | bc -l) )) && (( $(echo "$bBRO_pval > 0.005" | bc -l) ));
        then echo $OG >> DE.in.BGM_no.BRO.lst; echo $OG
        fi;

done < orthogroups-de_with_m.ref.tab

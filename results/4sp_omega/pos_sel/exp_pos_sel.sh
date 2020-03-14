while read line;

do OG=$(echo $line | awk '{print $1}');

model_BAT=$(echo $line | awk '{print $2}');
model_BRO=$(echo $line | awk '{print $4}');
model_BGM=$(echo $line | awk '{print $6}');

	if [[ $model_BRO == alternative ]] && [[ $model_BAT == alternative ]] && [[ $model_BGM == general ]];
	then echo $OG >> automictic_pos_sel_genes.lst;
	fi;

        if [[ $model_BRO == alternative ]] && [[ $model_BAT == general ]] && [[ $model_BGM == general ]];
        then echo $OG >> BRO_pos_sel_genes.lst;
        fi;

        if [[ $model_BRO == general ]] && [[ $model_BAT == alternative ]] && [[ $model_BGM == general ]];
        then echo $OG >> BAT_pos_sel_genes.lst;
        fi;

done < 4sp_pos_sel.tab

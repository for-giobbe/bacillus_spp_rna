echo "OG BAT_bs_omega_model BAT_bs_omega_sites BRO_bs_omega_model BRO_bs_omega_sites BGM_bs_omega_model BGM_bs_omega_sites" > 4sp_omega_pos_sel.tab;

for i in 4_sp/*aln;

do

#	unset $og $BAT_model $BAT_sites $BRO_model $BRO_sites $BGM_model $BGM_sites

	og=$(echo $i | awk -F "/" '{print $2}' | awk -F "." '{print $1}');

        if ls 4_sp_analyze_BAT_tot/$og.ref.mafft.n_replicate_?_model_alternative.out || ls 4_sp_analyze_BAT_tot/$og.ref.mafft.n_replicate_?_model_general.out;
        then
	BAT_model=$(echo 4_sp_analyze_BAT_tot/$og*.out | awk -F "_" '{print $NF}' | awk -F "." '{print $1}' );
	BAT_sites=$(awk 'x==1 {print $0} /Bayes Empirical Bayes/ {x=1}' 4_sp_analyze_BAT_tot/$og*.out | grep -c "*");
	else
    	BAT_model=NA;
        BAT_sites=NA;
        fi;


	if ls 4_sp_analyze_BRO_tot/$og.ref.mafft.n_replicate_?_model_alternative.out || ls 4_sp_analyze_BRO_tot/$og.ref.mafft.n_replicate_?_model_general.out;
        then
	BRO_model=$(echo 4_sp_analyze_BRO_tot/$og*.out | awk -F "_" '{print $NF}' | awk -F "." '{print $1}' );
	BRO_sites=$(awk 'x==1 {print $0} /Bayes Empirical Bayes/ {x=1}' 4_sp_analyze_BRO_tot/$og*.out | grep -c "*");
	else
	BRO_model=NA;
        BRO_sites=NA;
        fi;


	if ls 4_sp_analyze_BGM_tot/$og.ref.mafft.n_replicate_?_model_alternative.out || ls 4_sp_analyze_BGM_tot/$og.ref.mafft.n_replicate_?_model_general.out;
	then
	BGM_model=$(echo 4_sp_analyze_BGM_tot/$og*.out | awk -F "_" '{print $NF}' | awk -F "." '{print $1}' );
	BGM_sites=$(awk 'x==1 {print $0} /Bayes Empirical Bayes/ {x=1}' 4_sp_analyze_BGM_tot/$og*.out |  grep -c "*");
	else
	BGM_model=NA;
	BGM_sites=NA;
	fi;

	echo "$og $BAT_model $BAT_sites $BRO_model $BRO_sites $BGM_model $BGM_sites" >> 4sp_omega_pos_sel.tab;

done

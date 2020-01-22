
	for i in *aln;

		do OG=$(echo $i | awk -F "." '{print $1}'); printf $OG;

		for j in $(grep ">T" $i);

			do transcript=$(echo $j | awk -F "_i" '{print $1}' | tr -d ">");

			species=$(echo $j | awk -F "_" '{print $NF}');

			t=$(grep $transcript $species".result" | awk '{print $1" "$6" "$NF}');

			printf " $species  $t" >> test.tmp;

		done ;

	printf " $OG \n" >> test.tmp; echo -e "";

	done;



	while read line;

		do echo $line | awk '/BAT TR*/ && /BRO TR*/ && /BGM TR*/';

	done < test.tmp > orthogroups-de_wout_m.tab;



	while read line;

		do transcript=$(echo $line | awk '{print $10}');

		missing_line=$(grep $transcript BGM_m_only.results | awk '{print $1" "$6" "$NF}');

		echo BGM_M $missing_line $line >> orthogroups-de_with_m.tmp;

	done < orthogroups-de_wout_m.tab

	while read line; 

		do echo $line | awk '/BGM_M TR*/ && /BAT TR*/ && /BRO TR*/ && /BGM TR*/'; 

	done < orthogroups-de_with_m.tmp > orthogroups-de_with_m.tab

        rm *tmp


	for i in *aln;

		do OG=$(echo $i | awk -F "." '{print $1}'); printf $OG;

		for j in $(grep ">" $i);

			do transcript=$(echo $j | awk -F "_i" '{print $1}' | tr -d ">");

			species=$(echo $j | awk -F ".p[0-9]" '{print $2}' | tr -d "_1");

			t=$(grep $transcript $species".result" | awk '{print $1" "$6" "$NF}');

			printf " $species  $t" >> test.tmp;

		done ;

	printf " $OG \n" >> test.tmp; echo -e "";

	done;



	while read line; 

		do echo $line | awk '/BAT TR*/ && /BRO TR*/ && /BGM TR*/'; 

	done < test.tmp > orthogroups-de.tab;

	rm *test.tmp

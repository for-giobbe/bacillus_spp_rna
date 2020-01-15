for i in $(awk -F "\t" '{print $1}' $1 | sort -u);

        do

	a=$(grep $i $1 | awk -F "\t" '{print $(NF-1)" "$NF}' | sort -g | head -500 | awk '{print $NF}');

        echo $a | sed "s/ /\n/g" | sed "s/;/\n/g" > $1"tmp.taxon.lst";

 cat $1"tmp.taxon.lst"

        taxonkit lineage $1"tmp.taxon.lst" > $1"tmp.lineage.lst" ;

 cat $1"tmp.lineage.lst"

        tot=$(wc -l $1"tmp.lineage.lst" | awk '{print $1}');

        if [[ $(grep -c Panarthropoda $1"tmp.lineage.lst") -ge tot ]] #&& [[ $(grep -c Panarthropoda $1"tmp.lineage.lst") -ge 1 ]];

                then echo $i >> $3;

                else echo $i >> $2;

        fi;

done

echo "contaminants" >> $2
echo "insects" >> $3
#




head_first=$(head -1 $1);

head_second=$(head -1 $2);

echo $head_first $head_second > $3.tab;

tail -n +2 $1 > tab.tmp

while read line;

        do

	og=$(echo $line | awk '{print $1}');

        newline=$(grep $og $2);

        echo "$line $newline" >> $3.tab;

        done < tab.tmp

rm tab.tmp

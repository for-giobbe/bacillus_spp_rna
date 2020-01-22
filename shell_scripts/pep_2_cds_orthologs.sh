declare -a arr=("$@");

for file in *.fa;

        do

	og=$(echo $file | awk -F "\." '{print $1}');

        for sp in "${arr[@]}" ;

        	do seq_name=$(grep $sp $file | awk -F "_" 'sub(FS $NF,x)' | sed "s/\.p.[0-9]*//");

                echo $seq_name"_"$sp >> $og".cds.fa";

                grep -A 1 $seq_name $sp".final.raw.transdecoder.cds.ol" | tail -1 >> $og".cds.fa";

                done;

        done

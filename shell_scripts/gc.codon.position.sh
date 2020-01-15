
while getopts ":c:i:o:h" o; do
    case "${o}" in

        c) c=${OPTARG}
            ;;
	i) i=${OPTARG}
            ;;
	h) echo "
                        This script will calculate GC contant at 1+2 and 3rd codon positions for each contig and then sort them in two bin on the basis of GC3.
			
			List of non-optional arguments:
			-i fasta-formatted input file.
			-c GC3 cutoff (a number between 1 and 99).
			-h help page.

			This script requires infoseq and seqkit, it has been teste with conda installations.
"
               exit
           ;;
       \?) echo "WARNING! -$OPTARG isn't a valid option"
           exit
          ;;
       :) echo "WARNING! missing -$OPTARG value"
          exit
          ;;
       esac
 done

if [ -z "$c" ] || [ -z "$i" ]
then
 echo " WARNING! non-optional argument/s is missing "
exit
fi

#####################################################################################################################


name=$(echo $i | awk -F "\." '{ $NF="";print}' | sed "s/ /\./g")

echo -e "gene 12 03 tot" > $name"gc"

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < $i > oneliner.gc.tmp

while read line; do if echo $line | grep -qv ">"; then echo $line | sed "s/\(..\)./\1/g" >> 1+2.position.gc.tmp; else echo $line >> 1+2.position.gc.tmp; fi; done < oneliner.gc.tmp

while read line; do if echo $line | grep -qv ">"; then echo $line | sed 's/.\{2\}\(.\)/\1/g' >> 3rd.position.gc.tmp; else echo $line >> 3rd.position.gc.tmp; fi; done < oneliner.gc.tmp

infoseq -auto -pgc 1+2.position.gc.tmp | awk -F " " '{print $1" "$7}' | awk -F ":" '{print $4}' | tail -n +2 > 1+2.gc.tmp

infoseq -auto -pgc 3rd.position.gc.tmp | awk -F " " '{print $1" "$7}' | awk -F ":" '{print $4}' | tail -n +2 > 3rd.gc.tmp

infoseq -auto -pgc oneliner.gc.tmp | awk -F " " '{print $1" "$7}' | awk -F ":" '{print $4}' | tail -n +2 > tot.gc.tmp

paste 1+2.gc.tmp 3rd.gc.tmp tot.gc.tmp | awk '{print $1" "$2" "$4" "$6}' > $name"gc"

rm *gc.tmp

while read line; 

	do gc=$(echo $line | awk '{print $3}' | awk -F "\." '{print $1}'); 

		if [[ $gc -ge $c ]]; then echo $line | awk '{print $1}' >> $name"GC3>"$c".lst"; 

		else echo $line | awk '{print $1}' >> $name"GC3<"$c".lst"	

		fi; 

done < $name"gc"

seqkit grep -f $name"GC3>"$c".lst" $i > $name"GC3>"$c".fasta"

seqkit grep -f $name"GC3<"$c".lst" $i > $name"GC3<"$c".fasta"


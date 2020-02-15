#first argument is vcf
#second argument is number of readgroups

n=$(($2+9))

#hom=$(for i in $(seq 10 "$n"); do awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "0/0"; done);
#het=$(for i in $(seq 10 "$n"); do awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "0/1\|1/0"; done);
#uid=$(for i in $(seq 10 "$n"); do awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "\./\."; done);


echo "hom het NA perc tot"

for i in $(seq 10 "$n"); do

	hom=$(awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "0/0"); 
	het=$(awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "0/1\|1/0");
	uid=$( awk -vcol="$i" -v OFS="\t" '$1 !~ "^#" {print $col}' $1 | grep -c "\./\.");
	tot=$((het+hom+uid))
	per=$(awk "BEGIN { pc=100*${hom}/${tot}; i=int(pc); print (pc-i<0.5)?i:i+1 }")
	echo "$hom $het $uid $per $tot"
done
#echo -n " $(seq 1 "$2") \n" | tr ' ' \\t
echo -e hom sites $hom | tr ' ' \\t > $1.variants
echo -e het sites $het | tr ' ' \\t >>$1.variants
echo -e unident.  $uid | tr ' ' \\t >>$1.variants

# arguments can be euther BGM_S BRO_F BAT_F BGM_M BGM_F

echo "transcript_$1 len_$1 n_var_$1" > $1"_var.tab";

sp=$(echo $1 | awk -F "_" '{print $1}')

echo $sp

for i in $(grep "^TRINITY" $1".recode.vcf" | awk '{print $1}' | sort -u);

	do

	n_var=$(grep -c $i $1".recode.vcf");

	len=$(grep -A 1 $i $sp".def.nuc.fa" | tail -1 | wc -c)

	echo -e "$i $len $(($n_var -1 ))" >> $1"_var.tab";

	done

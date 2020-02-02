for i in *aln;

do

	name=$(echo $i | awk -F "." '{print $1}');

	Gblocks $i -t c -b5 h;

	rm *htm;

	mv $i Sname.ref.mafft.n.gb.aln;

done

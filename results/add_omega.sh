while read line; 
 do OG=$(echo $line | awk '{print $NF}'); echo -n $line >> test.tmp; 
 for i in BAT BRO BGM BGM_BAT; do newline=$(grep -w $OG 4sp_omega/branch.$i.*dNdS.summary); 
 echo -n " $newline" >> test.tmp; done; echo -e "" >> result.tab; done < results.tab
  sed -i 's/ \t / /g' result.tab

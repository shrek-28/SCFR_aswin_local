#Filter SCFRs that are part of coding region

#Filter with different length cut offs (16.0167 mins)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for win in 100 500 1000 2500 5000 7500 10000
do
echo ">"$win
mkdir -p /media/aswin/SCFR/SCFR-main/SCFR_lists/"$win"
  for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo " -"$species
  cds=$(find genes/$species/ -name "GCF_*.bed")
  #Filter SCFRs longer than the window
  awk -v a="$win" '{if($3-$2>=a) print$1,$2,$3,$4}' OFS="\t" SCFR_all/${species}_SCFR_all.out > SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out"
  #Get SCFRs that don't overlap with coding genes
  bedtools intersect -v -a SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" -b $cds | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed"
  ) &
  done
wait
unset species cds
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e


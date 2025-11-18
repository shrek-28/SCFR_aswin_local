#Filter SCFRs that are part of coding region

#mv SCFR_lists github/

mkdir -p /media/aswin/SCFR/SCFR-main/SCFR_lists/5kb
cd /media/aswin/SCFR/SCFR-main

start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo $species
#mkdir -p /media/aswin/SCFR/SCFR-main/SCFR_lists/$species
cds=$(find genes/$species/ -name "GCF_*.bed")
#Filter SCFRs longer than 5kb
awk '{if($3-$2>=5000) print$1,$2,$3,$4}' OFS="\t" SCFR_all/${species}_SCFR_all.out > SCFR_lists/5kb/${species}"_SCFR_atleast_5kb.out"
#Get SCFRs that don't overlap with coding genes
bedtools intersect -v -a SCFR_lists/5kb/${species}"_SCFR_atleast_5kb.out" -b $cds > SCFR_lists/5kb/$species"_scfr_atleast_5kb_in_non_coding.bed"
#Get SCFRs that don't overlap with coding genes
bedtools intersect -v -a SCFR_all/${species}_SCFR_all.out -b $cds > SCFR_lists/5kb/$species"_scfr_in_non_coding.bed"
#Filter SCFRs longer than 5kb
awk '{if($3-$2>=5000) print$1,$2,$3,$4,$3-$2}' OFS="\t" SCFR_lists/5kb/$species"_scfr_in_non_coding.bed" > SCFR_lists/5kb/$species"_scfr_in_non_coding_atleast_5kb.bed"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

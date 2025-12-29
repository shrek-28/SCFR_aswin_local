
cd /media/aswin/SCFR/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
mkdir exon_shadow/"$species"
cd exon_shadow/"$species"
#get zero based scfr bed file
time awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/"$species"_SCFR_all.out > "$species"_scfr_all.bed
#Convert gtf to zero based cds bed (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
gtf=$(find /media/aswin/SCFR/SCFR-main/genes/"$species" -name "GCF*.gtf")
time python3 /media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/gtf_to_cds_with_transcript_exon_metadata.py -i $gtf -s "$species" -o "$species"_coding_exons.bed
#get SCFR exon overlaps (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time bedtools intersect -a "$species"_scfr_all.bed -b "$species"_coding_exons.bed -wo -s > "$species"_scfr_cds_all_overlaps.bed
#get scfrs containing exons with same frame as exons & calculate their upstream & downstream shadow
time awk '{c4=$4; gsub("-","",c4); if(c4==$13 && $2<=$8 && $3>=$9) {if($6=="+") print$0,$8-$2,$3-$9; else if($6=="-") print$0,$3-$9,$8-$2}}' "$species"_scfr_cds_all_overlaps.bed > "$species"_scfr_containing_cds.bed
#filter & classiffy exon-scfr overlaps & calculate exon shadow
#create files
echo "chr start end frame filler strand chr exon_start exon_end gene transcript exon_strand exon_frame exon_number gene_tx_count exon_tx_count exon_order exon_sharing exon_splicing" | tr " " "\t" > "$species"_single_exon.tsv
echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap" | tr " " "\t" > "$species"_multi_exon.tsv
#Loop over each unique SCFRs
time while read scfr
do
grep "$scfr" "$species"_scfr_containing_cds.bed > scfr_temp.bed
#strand
strand=$(awk -F "\t" '{print$6}' scfr_temp.bed | sort -u)
#Total number of exons
tnoe=$(wc -l < scfr_temp.bed)
#Number of overlapping exons
mc=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" scfr_temp.bed | bedtools sort -i - | bedtools merge -i - -s -c 5 -o count | wc -l)
#Classify exon shadow
if [[ $tnoe == 1 ]] || [[ $mc == 1 ]]
then
#save unique single exon containing SCFRs (single exon per SCFR)
#filter identical exons, keep only one
awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' scfr_temp.bed >> "$species"_single_exon.tsv
elif [[ $mc > 1 ]]
then
#save multi-exon containing SCFRs (multiple exons per SCFR or single transcript per SCFR)
for ut in $(awk '{print$11}' scfr_temp.bed | sort -u)
do
awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8n -k9,9n > temp_multi_sorted
#plus strand
if [[ $strand == "+" ]]
then
e1=$(head -1 temp_multi_sorted | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}')
e2=$(tail -1 temp_multi_sorted| awk '{print$9,$14,$16,$17,$18,$19,$20,$22}')
ec=$(wc -l < temp_multi_sorted)
paste <(echo "$e1") <(echo "$e2") | awk -v ec="$ec" '{print$1,$2,$3,$4,$5,$6,$7,$8,$21,$9,$10,$11,$12,ec,$20,$28,$13,$22,$14,$15,$23,$16,$24,$17,$25,$18,$26,$19,$27}'
#minus strand
elif [[ $strand == "-" ]]
then
e1=$(head -1 temp_multi_sorted | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$22}')
e2=$(tail -1 temp_multi_sorted | awk '{print$9,$14,$16,$17,$18,$19,$20,$21}')
ec=$(wc -l < temp_multi_sorted)
paste <(echo "$e1") <(echo "$e2") | awk -v ec="$ec" '{print$1,$2,$3,$4,$5,$6,$7,$8,$21,$9,$10,$11,$12,ec,$28,$20,$22,$13,$14,$23,$15,$24,$16,$25,$17,$26,$18,$27,$19}'
else :
fi
rm temp_multi_sorted
unset e1 e2 ec
#Filter transcripts with identical exons, keep only one
done | awk '!seen[$7,$8,$9,$10,$12,$13,$15,$16,$22,$23,$26,$27,$28,$29]++' | tr " " "\t" >> "$species"_multi_exon.tsv
unset ut
else :
fi
#delete temp files & unset variables
rm scfr_temp.bed
unset tnoe mc
done < <(awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" "$species"_scfr_containing_cds.bed | sort -u)
unset scfr
cd /media/aswin/SCFR/SCFR-main
done


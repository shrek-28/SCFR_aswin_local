

bedtools intersect -a scfr.bed -b cds.bed -wo -s > scfr_cds_overlaps.bed

while read scfr 
do
f=$(echo $scfr | awk '{print$5}'
grep "$scfr" scfr_cds_overlaps.bed | awk '$5==$11 && $2<=$8 && $3>=$9 {print$0,$8-$2,$3-$9}'
done < scfr.bed


#Refomrat 
#1m58.476s
time awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/human_SCFR_all.out > human_scfr_all.bed
awk '{if($6~"-") print$1,$2,$3,$4,$6,"-",$5,$7,$8,$9,$10,$11,$12; else print$1,$2,$3,$4,$6,"+",$5,$7,$8,$9,$10,$11,$12}' OFS="\t" ../human_coding_exons.bed > human_coding_exons_reformatted.bed
#Get overlaps ()
time bedtools intersect -a human_scfr_all.bed -b human_coding_exons_reformatted.bed -wo -s > human_scfr_cds_all_overlaps.bed

#Convert gtf to cds bed (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time python3 gtf_to_cds_with_transcript_exon_metadata_2.py -i /media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf -o human_coding_exons.bed
#Get overlaps (2m14.940s) (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
time bedtools intersect -a human_scfr_all.bed -b human_coding_exons.bed -wo -s > human_scfr_cds_all_overlaps.bed

while read scfr 
do
f=$(echo $scfr | awk '{print$4}' | tr -d "-")
#restrict exon overlaps to within scfr or at scfr boundaries & have same frame
grep "$scfr" human_scfr_cds_all_overlaps.bed | awk '{gsub(/-/,"",$4)}1' OFS="\t" | awk '$4==$13 && $2<=$8 && $3>=$9 {print$0}' > temp.bed
#Total number of exons
tnoe=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" temp.bed | sort -u | wc -l)
#Overlapping exons
mc=$(bedtools sort -i temp.bed | bedtools merge -i - | wc -l)
if [[ $tnoe == "1" ]] || [[ $mc == "1" ]]
then
cat temp.bed | awk '{if($6~"-") print$0,$3-$9,$8-$2; else print$0,$8-$2,$3-$9}' OFS="\t" | awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' >> single_exon.tsv

c=$(bedtools sort -i temp.bed | wc -l)




grep "$scfr" human_scfr_cds_all_overlaps.bed | awk '$4==$11 && $2<=$8 && $3>=$9 {print$0,$8-$2,$3-$9}'
done < <(awk '{print$1,$2,$3,$4,$5,$6}' human_scfr_cds_all_overlaps.bed | sort -u | tr " " "\t")


scfr=`grep KLHL17 human_scfr_cds_all_overlaps.bed | grep  392806 | grep  393139 | grep XM_054336255.1 | grep last`
scfr=`grep BRCA1 human_scfr_cds_all_overlaps.bed | grep  43926117 | grep 43926252`
scfr=`grep TTC23L human_scfr_cds_all_overlaps.bed | grep  35087777 | grep 35088032 | awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" | sort -u`

scfr=`grep KLHL17 human_scfr_cds_all_overlaps.bed | grep 393139 | grep 392806 | grep XM_054336255.1 | grep 392897`



awk '!seen[$7,$8,$9,$10,$12,$13,$17,$19,$20,$21,$22]++' 

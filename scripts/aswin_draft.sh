#############################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
#############################################################################################################################################################################################################################################################################################################


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#BLAST SCFR against nr database

cd /media/aswin/SCFR/SCFR-main/gene_deserts/SCFR_overlap_gene_deserts/ncbi_nr_search
time blastp -query test.fa -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/v5_nr_blastdb/nr -out results.blastp.out -evalue 0.001 -max_target_seqs 100 -outfmt 6 -num_threads 32

cd /media/aswin/SCFR/SCFR-main/
start_time=$(date +%s)
#Don't use gorilla here
for species in human chimpanzee bonobo gibbon borangutan sorangutan
do
echo ">"$species
cd gene_deserts/SCFR_overlap_gene_deserts/
mkdir $species/SCFR_fasta/v4
for path in $(find $species/SCFR_fasta -name "*_overlapping_scfrs_canonical_orf.fa")
do
can=$(echo $path)
noncan=$(echo $path | sed 's/_canonical_/_non_canonical_/g')
scfrcan=$(echo $path | awk -F "/" '{print$NF}')
echo " - "$scfrcan
#Run blast on canonical ORFs
	#time blastp -query $can -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.001 -max_target_seqs 20 --max_hsps 3 -outfmt 11 -num_threads 32 -out $species/SCFR_fasta/$scfrcan".outfmt11"
	#blast_formatter -archive $species/SCFR_fasta/$scfrcan".outfmt11" -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
		# | sed '1i Sub_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > $species/SCFR_fasta/$scfrcan".outfmt6"
time blastp -query $can -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.0001 -max_target_seqs 20 -max_hsps 3 -qcov_hsp_perc 70 -num_threads 32 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
 | sed '1i Subject_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > $species/SCFR_fasta/v4/$scfrcan".outfmt6"

#Run blast on non-canonical ORFs
fnoncan=$(echo $noncan | sed 's/_non_canonical_orf.fa/_non_canonical_filtered_orf.fa/g')
python3 /media/aswin/SCFR/SCFR-main/my_scripts/subtract_fasta.py $can $noncan > $fnoncan
scfrnoncan=$(echo $fnoncan | awk -F "/" '{print$NF}')
	#time blastp -query $fnoncan -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.001 -max_target_seqs 20 -max_hsps 3 -outfmt 11 -num_threads 32 -out $species/SCFR_fasta/$scfrnoncan".outfmt11"
	#blast_formatter -archive $species/SCFR_fasta/$scfrnoncan".outfmt11" -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
		#| sed '1i Sub_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > $species/SCFR_fasta/$scfrnoncan".outfmt6"
time blastp -query $fnoncan -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.0001 -max_target_seqs 20 -max_hsps 3 -qcov_hsp_perc 70 -num_threads 32 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
 | sed '1i Sub_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > $species/SCFR_fasta/v4/$scfrnoncan".outfmt6"
unset can noncan scfrcan fnoncan scfrnoncan
done
cd /media/aswin/SCFR/SCFR-main/
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e



time blastp -query test.fa -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.001 -max_target_seqs 100 -outfmt 11 -num_threads 32 -out test.outfmt11
blast_formatter -archive test.outfmt11 -outfmt 3 -line_length 280 -out test.outfmt3
blast_formatter -archive test.outfmt11 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
| sed '1i Sub_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > test.outfmt6

time blastp -query $fnoncan -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4/nr -evalue 0.001 -max_target_seqs 20 -max_hsps 3 -num_threads 32 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
 | sed '1i Sub_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | tr " " "_" | column -t > $species/SCFR_fasta/$scfrnoncan".outfmt6"


# Convert your protein database (nr) into a DIAMOND format
diamond makedb --in blast_nr_v4/nr.fasta -d nr_diamond_db

#Get fasta from pre-existing databases
cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/blast_nr_v4
blastdbcmd -db nr -entry all -outfmt %f -out nr.fasta
diamond makedb --in nr.fasta --db nr_diamond_db

diamond blastp \
-q input.fa \
-d nr_diamond_db \
-o SCFR_protein_hits.tsv \
-e 0.001 \
--max-target-seqs 20 \
-p 32 \
--outfmt 6 \
--more-sensitive


diamond blastp \
-q input.fa \
-d nr_diamond_db \
-o SCFR_protein_hits_sensitive.tsv \
-e 0.001 \
--max-target-seqs 100 \
--max-hsps 3 \
-p 32 \
--more-sensitive \
--matrix BLOSUM45 \
--outfmt 6
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#

  mkdir /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
  cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
  
  cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
  find /media/aswin/SCFR/SCFR-main/genes/$species/ -name "*merged*" | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot/'
  report=$(readlink -f /media/aswin/SCFR/SCFR-main/genome_reports/* | grep $species)
  awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1'  $report | sed 's/ /_/g' | awk 'NR>1{print$9,$12}' | sed '1i accession names' > $species"_chromosome_names"
  cp /media/aswin/SCFR/SCFR-main/genome_sizes/$species".genome" .


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#View SCFR length filtered summary in window-wise

for i in $(awk -F "," '{print$2}' window_wise_all_species_scfr_coding_stats.csv | grep -v Window | sort -u)
do
j=$(for s in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk -F "," -v a="$i" '$2==a {print$1,$3}' filtered_window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sort -k1,1n | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' | \
awk 'NR==1{print; next} {for(i=2;i<=NF;i++) {if($i>=100000) {$i=sprintf("%.2fm",$i/1000000)} else if($i>=1000) {$i=sprintf("%.2fk",$i/1000)} else {$i=sprintf("%.0f",$i)}} print}' OFS="\t" | column -t > length_threshold_total_scfr_count_summary

for i in $(awk -F "," '{print$2}' window_wise_all_species_scfr_coding_stats.csv | grep -v Window | sort -u)
do
j=$(for s in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk -F "," -v a="$i" '$2==a {print$1,$7}' filtered_window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sort -k1,1n | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' | column -t > length_threshold_Percent_genome_covered_by_SCFR_summary

for i in $(awk -F "," '{print$2}' window_wise_all_species_scfr_coding_stats.csv | grep -v Window | sort -u)
do
j=$(for s in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk -F "," -v a="$i" '$2==a {print$1,$6}' filtered_window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sort -k1,1n | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' | column -t > length_threshold_Percent_unfiltered_by_filtered_SCFR_summary

for i in $(awk -F "," '{print$2}' window_wise_all_species_scfr_coding_stats.csv | grep -v Window | sort -u)
do
j=$(for s in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk -F "," -v a="$i" '$2==a {print$1,$8}' filtered_window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sort -k1,1n | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' | column -t > length_threshold_Percent_SCFR_by_coding_summary

for i in $(awk -F "," '{print$2}' window_wise_all_species_scfr_coding_stats.csv | grep -v Window | sort -u)
do
j=$(for s in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk -F "," -v a="$i" '$2==a {print$1,$13}' window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sort -k1,1n | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' \
| awk 'NR==1{print; next} {for(i=2;i<=NF;i++) {if($i>=100000) {$i=sprintf("%.2fm",$i/1000000)} else if($i>=1000) {$i=sprintf("%.2fk",$i/1000)} else {$i=sprintf("%.0f",$i)}} print}' OFS="\t" | column -t > length_threshold_SCFR_cds_overlap_summary

cat length_threshold_percent_genome_covered_summary | sed -n 4p | awk '!($1="")' | tr " " "\n" | awk NF | sort -r

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#View higher strand asymmetry
sed 's/ /_/g' genome_reports/GCA_029289425.3_bonobo.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '{print$9,$3,$4,$5,$11}' | sort -k5,5n | colnum.sh 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#gene deserts scfr overlaps

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
echo ">"$species
mkdir gene_deserts/SCFR_overlap_gene_deserts/$species
for len in 0 100 500 1000 2500 5000 7500 10000
do
(
/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a SCFR_lists/$len/$species"_SCFR_atleast_"$len".out" -b gene_deserts/fishers_test/$species/$species"_only_intergenic_gene_deserts.bed" -wo > gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_"$len"_only_intergenic_gene_deserts_overlaps.out"
) &
done
wait
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
echo ">"$species
for len in 0 100 500 1000 2500 5000 7500 10000
do
(
/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a SCFR_lists/$len/$species"_SCFR_atleast_"$len".out" -b gene_deserts/fishers_test/$species/$species"_intronic_intergenic_gene_deserts.bed" -wo > gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_"$len"_intronic_intergenic_gene_deserts_scfr_overlaps.out"
) &
done
wait
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

/media/aswin/SCFR/SCFR-main/PCA/human/5000/with_coding_region/NC_060928.1.fasta

#DRAFT SCRIPTS
####################################################################################################################################################################################################################################################################################################################

#Download large genome files
datasets download genome accession $genome --include genome,gtf,seq-report --dehydrated --filename $genome.zip
unzip $genome.zip -d $genome
time datasets rehydrate --directory $genome

#Fourier ananlysis
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna.gz
gzip -d GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna.gz
python3 split_fasta_by_chunks.py GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna .
time python3 fft_motif_report_grouped.py GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna -o .
python split_fasta_by_chunks.py GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna CDS200
python fft_motif_analysis.py GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna -o output_folder

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis
cp -r ../PCA/* .
find . -name "*.tsv" | xargs rm
find . -name "*.pdf" | grep -v fft | xargs rm

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
for win in 1000 2500 5000 7500 10000
cd $species/5000/without_coding_region/
for scfr in $(ls *.fasta)
do
echo " -"$scfr
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o "output_"$scfr -t 32 $scfr
done
unset scfr
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

##

grep -v "^#" GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf | awk '{for (i=1;i<=NF;i++) if($i ~/gene_biotype/) print $(i+1)}' | sort | uniq -c > annotation_types

##
#Get full genomes for cds extraction (3m35.115s)
mkdir /media/aswin/SCFR/SCFR-main/genomes
cd /media/aswin/SCFR/SCFR-main
time while read i
do
genome=$(echo $i | awk '{print$1}')
species=$(echo $i | awk '{print$2}')
echo "-" 
unzip -o $genome".zip" -d genomes
mkdir genomes/$species
mv genomes/ncbi_dataset/data/$genome/* genomes/$species/
samtools faidx genomes/$species/GCA*.fna
unset genome species
done < QC/genome_accessions
rm -r genomes/ncbi_dataset genomes/md5sum.txt genomes/README.md

#Extract cds from genomes based on gtf
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
while read i
do
(
species=$(echo $i | awk '{print$2}')
genome=$(find genomes/$species -name "GCA*.fna")
gtf=$(readlink -f genes/$species/GC*.gtf)
echo -e ">"$species "\n -genome: "$genome"\n -GTF: "$gtf
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/$species
gffread $gtf -g $genome -x /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/$species/${species}"_cds.fa"
) &
done < QC/genome_accessions
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Get cds sequences from gtf
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo ">"species
gtf=$(readlink -f genes/$species/GC*.gtf)
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/$species
cat chrs/$species/*.fasta > chrs/$species/$species"_genome.fa"
samtools faidx chrs/$species/$species"_genome.fa"
time gffread $gtf -g chrs/$species/$species"_genome.fa" -x /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/$species/${species}"_cds.fa"
rm chrs/$species/$species"_genome.fa"
unset chr gtf
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

####################################################################
#get length stats of all scfrs

time for scfr in $(ls | grep "_SCFR_all.out")
do
echo ">"$scfr
python3 all_scfr_stats.py $scfr > "length_stats_"$scfr
done

####################################################################
#Window wise all species scfr cds stats

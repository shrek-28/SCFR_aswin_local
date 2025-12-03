#############################################################################################################################################################################################################################################################################################################



mkdir /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot

cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
find /media/aswin/SCFR/SCFR-main/genes/$species/ -name "*merged*" | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot/'
report=$(readlink -f /media/aswin/SCFR/SCFR-main/genome_reports/* | grep $species)
awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1'  $report | sed 's/ /_/g' | awk 'NR>1{print$9,$12}' | sed '1i accession names' > $species"_chromosome_names"
cp /media/aswin/SCFR/SCFR-main/genome_sizes/$species".genome" .


##########
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

#View higher strand asymmetry
sed 's/ /_/g' genome_reports/GCA_029289425.3_bonobo.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '{print$9,$3,$4,$5,$11}' | sort -k5,5n | colnum.sh 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

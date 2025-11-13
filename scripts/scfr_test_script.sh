####################################################################################################################################################################################################################################################################################################################
#
####################################################################################################################################################################################################################################################################################################################

####################################################################################################################################################################################################################################################################################################################

#Doownload github repository
cd /media/aswin/SCFR
wget https://github.com/ceglab/SCFR/archive/refs/heads/main.zip
unzip main.zip

#Create folders
cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  echo $species
  mkdir -p chrs/"$species"
  mkdir -p SCFR/"$species"
  mkdir -p SCFR/gaps/"$species"
  mkdir -p genes/"$species"
  mkdir -p GC/"$species"
  mkdir -p SCFR_all
  mkdir -p scripts
done

####################################################################################################################################################################################################################################################################################################################
#Quality check

cd /media/aswin/SCFR/SCFR-main
mkdir -p QC/genome_metadata

#Get genome metadata using datasets (0m18.381s)
cd /media/aswin/SCFR/SCFR-main/QC
while read i
  do
  echo $i
  i1=$(echo $i | awk '{print$1}')
  i2=$(echo $i | awk '{print$2}')
  datasets summary genome accession "$i1" --report sequence --as-json-lines > genome_metadata/$i2"_"$i1"_genome_metadata.json"
  cat genome_metadata/$i2"_"$i1"_genome_metadata.json" | json2xml | xtract -pattern opt -def "-" -element assembly_accession assembly_unit assigned_molecule_location_type chr_name gc_count gc_percent genbank_accession length refseq_accession role sequence_name | tr " " "_" | sed '1i assembly_accession assembly_unit assigned_molecule_location_type chr_name gc_count gc_percent genbank_accession length refseq_accession role sequence_name' | column -t > genome_metadata/$i2"_"$i1"_genome_metadata.txt"
  unset i1 i2
done < genome_accessions

cd /media/aswin/SCFR/SCFR-main/QC/genome_metadata
find . -name "*_genome_metadata.txt" | xargs -n1 sh -c 'awk "{print\$9,\$8}" $0' | grep -v "refseq_accession" | sed '1i refseq_accession Length' > all_species_chr_length

####################################################################################################################################################################################################################################################################################################################
#Download all the seven primate genomes

cd /media/aswin/SCFR/SCFR-main
time while read species
  do
  genome=$(echo $species | awk '{print$1}')
  name=$(echo $species | awk '{print$2}')
  echo ">"$genome $name
  #Download genome & annotation files
  #time datasets download genome accession $genome --include genome,gtf,seq-report --filename $genome".zip"
  unzip -o $genome".zip"
  cd ncbi_dataset/data/$genome
  #Get accessions
  cat sequence_report.jsonl | json2xml | xtract -pattern opt -def "-"  -element genbankAccession refseqAccession | awk '$2!="-"' > $genome"_genbank_refseq_accession"
  #Add refseq accesion in genome
  time awk 'NR==FNR {map[$1]=$2; next}  /^>/ {hdr=$0; for(k in map) {if(hdr ~ ">"k"($|[[:space:]])") {sub(">"k, ">"map[k]" "k); break}}} {print}' $genome"_genbank_refseq_accession" $genome*".fna" > $genome".fna"
  #Get chromosome fasta files
  outdir="/media/aswin/SCFR/SCFR-main/chrs/$name"
  awk -v outdir="$outdir" '/^>/ { if (out) close(out) ; out = outdir "/" substr($1,2) ".fasta" } { print >> out }' "${genome}.fna"
  unset genome name
  cd /media/aswin/SCFR/SCFR-main
done < QC/genome_accessions

####################################################################################################################################################################################################################################################################################################################
#QC: Check length of fetched genome sequences

cd /media/aswin/SCFR/SCFR-main/chrs
time for species in human bonobo chimpanzee gorilla borangutan gibbon
     do
     cd $species
     for chr in $(ls *.fasta)
     do
     length=$(awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $chr | awk '{print$1,$NF}')
     echo $species $length
     unset length
     done
     unset chr
     cd /media/aswin/SCFR/SCFR-main/chrs
done | sed '1i Species chromosome length' | column -t > /media/aswin/SCFR/SCFR-main/QC/genome_metadata/fetched_length

cd /media/aswin/SCFR/SCFR-main/QC/genome_metadata
for chr in $(awk '{print$1}' all_species_chr_length | egrep -v "refseq_accession|-")
     do
     m1=$(awk -v a="$chr" '$1==a' all_species_chr_length)
     m2=$(awk -v a="$chr" '$2==a' fetched_length)
     echo $m1 $m2
     unset m1 m2
done | awk '{if($2==$5) print$0,"same"; else print$0,"different"}' |  sed '1i ncbi_accession ncbi_length species fetched_accession fetched_length Compare_length' | column -t > compare_chromosome_length

####################################################################################################################################################################################################################################################################################################################
#List all the SCFRs in the genomes of the 7 primate species (47.0667 min without sorangutan)

#Time taken human (75m14.862s), bonobo (62m53.702s), chimpanzee (54m44.004s), gorilla (56m4.399s), borangutan (79m51.724s), sorangutan (70m45.710s), gibbon (66m39.910s)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  echo $species
  for chr in `ls -1 chrs/$species/*.fasta|cut -f 3 -d '/'|sed 's/\.fasta//g'`
  do
  (
  echo $chr
  python3 scripts/find_stop_codon_free_regions_with_reverse_gap_report.py chrs/$species/"$chr".fasta --gaps SCFR/gaps/"$species"/gap_regions_"$chr".bed > SCFR/"$species"/"$chr".fasta.SCFRs.out
  ) &
  done
  wait
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################
#collate all gaps

cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  echo $species
  cat SCFR/gaps/"$species"/gap*.bed|sort -k1,1 -k2n,2 > SCFR_all/gaps_"$species".bed
done

#########################################################################################################################
#get summary of all SCFRs in the genomes of the 7 primate species (251m11.981s)

#37.5167 min
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  (
  echo $species
  cat SCFR/"$species"/*.fasta.SCFRs.out | sort -k1,1 -k2n,2 | bedtools intersect -v -a stdin -b SCFR_all/gaps_"$species".bed > SCFR_all/"$species"_SCFR_all.out
  ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#38.7 min Since running 6 species need huge memory run 3 at a time
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  (
  echo $species
  Rscript scripts/summarize_SCFR_bed_frames_all.R SCFR_all/"$species"_SCFR_all.out
  ) &
  ((++i % 3 == 0)) && wait  # run 3 at a time
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################################################################################
#get strand assymetry of SCFRs per chromosome in the genomes of the 7 primate species (25m12.250s)

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  (
  echo $species
  python3 scripts/quantify_scfr_asymmetries_by_chrom.py SCFR_all/"$species"_SCFR_all.out SCFR_all/"$species"_SCFR_asymmetries_out.csv
  ) &
  done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#get strand assymetry of SCFRs in sliding windows across the genomes of the 7 primate species (46m48.751s)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan gibbon
  do
  (
  echo $species
  python3 scripts/quantify_scfr_asymmetries_by_chrom_window.py SCFR_all/"$species"_SCFR_all.out --window-size 100000 --slide-size 50000 --output SCFR_all/"$species"_SCFR_asymmetries_out_win100000_slide50000.csv
  ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#This script is actually saved as version_3_plot_strand_asymmetry_sliding_extremes.R in the scripts folder. The older versions were very sensitive to noise.
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo $species
  Rscript scripts/version_3_plot_strand_asymmetry_sliding_extremes.R --input SCFR_all/"$species"_SCFR_asymmetries_out_win100000_slide50000.csv --pdf SCFR_all/"$species"_sliding_outliers.pdf --bed SCFR_all/"$species"_sliding_outlier_regions.bed --min_region_size 100000 --window_len 5 --min_hits 3
  ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################################################################################
#get GC content of the SCFRs

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo $species
  cat chrs/$species/*.fasta > tmp."$species".fasta
  bedtools nuc -fi tmp."$species".fasta -bed SCFR_all/"$species"_SCFR_all.out > SCFR_all/"$species"_SCFR_GC_all.out
  rm tmp."$species".fasta
  cat SCFR_all/"$species"_SCFR_GC_all.out|awk '$13>10000{print $0}'|grep -v "^#"|sed 's/:/\t/g'|cut -f 1-13|sort -k1,1 -k2n,2 > SCFR_all/"$species"_long_SCFRs.bed
  ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################################################################################
##Download all the seven primate genome annotation files

cd /media/aswin/SCFR/SCFR-main
cd genes/human
#Download the gene annotation file of T2T-CHM13v2.0 genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
cd ../../
cd genes/chimpanzee
#Download the NHGRI_mPanTro3-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/bonobo
#Download the NHGRI_mPanPan1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/gibbon
#Download the NHGRI_mSymSyn1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf.gz
cd ../../
cd genes/gorilla
#Download the NHGRI_mGorGor1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gtf.gz
cd ../../
cd genes/borangutan
#Download the NHGRI_mPonPyg2-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/sorangutan
#Download the NHGRI_mPonAbe1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gtf.gz
cd ../../

#########################################################################################################################

#Extract coding exons annotated in the GTF
cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
gunzip genes/"$species"/*.gtf.gz
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
#Extract coding exon coordinates from GTF
python scripts/gtf_to_bed.py genes/"$species"/"$gtf".gtf genes/"$species"/"$gtf".bed
done

#########################################################################################################################
#bin the SCFR counts by length and GC content

cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
cat "$species"_SCFR_GC_all.out|awk '{printf "%s %.2f\n", $13, $5}' |sed 's/ /\t/g'|awk '{
  if ($1 < 1000)
    bin1 = int($1 / 100) * 100;
  else
    bin1 = int($1 / 1000) * 1000;

  bin2 = sprintf("%.1f", int($2 * 10) / 10);

  count[bin1, bin2]++;
  range[bin1] = ($1 < 1000) ? 100 : 1000;
} END {
  for (key in count) {
    split(key, bins, SUBSEP);
    r = range[bins[1]];
    printf "%d-%d\t%s-%s\t%d\n", bins[1], bins[1]+r-1, bins[2], sprintf("%.1f", bins[2]+0.1), count[key];
  }
}' > "$species"_bins.out
done

#Plot a 2 dimensional histogram of length vs AT content and label the SCFR longer than 10 Kb that overlap coding exons
cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
bedtools intersect -a SCFR_all/"$species"_long_SCFRs.bed -b genes/"$species"/"$gtf".bed -wa -wb|cut -f 5,13,17|sort -u > SCFR_all/"$species"_genes_of_interest.txt
done
cp scripts/combined_2dhist.r SCFR_all
cd SCFR_all
Rscript combined_2dhist.r

#########################################################################################################################












####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
####################################################################################################################################################################################################################################################################################################################

#Download
time datasets download genome accession GCA_028878055.3 --include genome,gtf,seq-report --filename GCA_028878055.3.zip

unzip -o GCA_028878055.3.zip
cd ncbi_dataset/data/GCA_028878055.3

#Get accessions
cat sequence_report.jsonl | json2xml | xtract -pattern opt -def "-"  -element genbankAccession refseqAccession | awk '$2!="-"' > genbank_refseq_accession

#Add refseq accesion in genome
awk 'NR==FNR {map[$1]=$2; next}  /^>/ {hdr=$0; for(k in map) {if(hdr ~ ">"k"($|[[:space:]])") {sub(">"k, ">"map[k]" "k); break}}} {print}' genbank_refseq_accession test2.fa > test2_updated.fna

#Add refseq accesion in genome
time awk 'NR==FNR {map[$1]=$2; next} /^>/ {hdr=$0; for(k in map) {if(hdr~">"k"($|[[:space:]])") {sub(">"k, ">"map[k]" "k); break}}} {print}' genbank_refseq_accession test1.fa > tmp && mv tmp test1.fa

time awk -iinplace 'NR==FNR {map[$1]=$2; next} /^>/ {id=substr($1,2); if (id in map) $0=">"map[id]" "id substr($0,length($1)+1)} {print}' genbank_refseq_accession test2.fna> tmp && mv tmp test2.fna
time gawk -iinplace 'NR==FNR {map[$1]=$2; next} /^>/ {id=substr($1,2); if (id in map) $0=">"map[id]" "id substr($0,length($1)+1)} {print}' genbank_refseq_accession test3.fna > tmp && mv tmp test2.fna

awk 'NR==FNR {map[$1]=$2;next} /^>/{hdr=$0 for(k in map) {if(hdr~">"k" ($|[[:space:]])") {sub(">"k, ">"map[k]" "k) break }}} {print}' genbank_refseq_accession GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna \
> updated.fna

start_time=$(date +%s)
while read chr 
do
(
g=$(echo $chr | awk '{print$1}')
r=$(echo $chr | awk '{print$2}')
echo $r $g
sed "/$g/ s/^>/>$r /g" test2.fna -i
) &
done < genbank_refseq_accession
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e


#
start_time=$(date +%s)

awk 'NR==FNR {map[$1]=$2; next}
     /^>/ {
         id=substr($1,2)
         if (id in map)
             $0 = ">" map[id] " " id substr($0, length($1)+1)
     }
     {print}' genbank_refseq_accession test2.fna > tmp && mv tmp test2.fna

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Elapsed time: ${elapsed_time}s"

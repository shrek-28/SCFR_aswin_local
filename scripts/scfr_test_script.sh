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
find . -name "*_genome_metadata.txt" | xargs -n1 sh -c 'awk "{print\$9,\$8}" $0' | grep -v "refseq_accession" | sed '1i refseq_accession Length' | grep -v "^\-" > all_species_chr_length

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
  #If genome size is too big or the downloads gets failed, use dehydrate
    #datasets download genome accession $genome --include genome,gtf,seq-report --dehydrated --filename $genome.zip
    #unzip $genome.zip -od $genome
    #datasets rehydrate --directory $genome
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

#Delete mitochondrial genomes
cd /media/aswin/SCFR/SCFR-main
grep -vf <(awk '{print$1}' QC/genome_metadata/all_species_chr_length | grep -v "refseq_accession") <(find chrs -name "*.fasta") | xargs rm

####################################################################################################################################################################################################################################################################################################################
#QC: Check length of fetched genome sequences

cd /media/aswin/SCFR/SCFR-main/chrs
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
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
#List all the SCFRs in the genomes of the 7 primate species (51.6833 mins)

#Time taken human (75m14.862s), bonobo (62m53.702s), chimpanzee (54m44.004s), gorilla (56m4.399s), borangutan (79m51.724s), sorangutan (70m45.710s), gibbon (66m39.910s)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
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
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  echo $species
  cat SCFR/gaps/"$species"/gap*.bed|sort -k1,1 -k2n,2 > SCFR_all/gaps_"$species".bed
done

#########################################################################################################################
#get summary of all SCFRs in the genomes of the 7 primate species (251m11.981s)

#(35.4167 mins)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo $species
  cat SCFR/"$species"/*.fasta.SCFRs.out | sort -k1,1 -k2n,2 | bedtools intersect -v -a stdin -b SCFR_all/gaps_"$species".bed > SCFR_all/"$species"_SCFR_all.out
  ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Since running 6 species need huge memory run 3 at a time (1.67528 hours ~100.517 mins)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
max_jobs=3
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo $species
  Rscript scripts/summarize_SCFR_bed_frames_all.R SCFR_all/"$species"_SCFR_all.out
  ) &
  # If max_jobs reached, wait for ONE job to finish (not all)
  if (( $(jobs -r | wc -l) >= max_jobs )); then
  wait -n
  fi
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e
unset max_jobs

#########################################################################################################################
#get strand assymetry of SCFRs per chromosome in the genomes of the 7 primate species (25m12.250s)

#11.3 mins
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
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
#9.45 mins
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
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
#6 secs
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

##################################################################################################################################################################################################################################################
#get GC content of the SCFRs

#48.1667
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

#6.18333 mins
cd /media/aswin/SCFR/SCFR-main/genes
start_time=$(date +%s)
declare -A urls=(
  [human]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
  [chimpanzee]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gtf.gz"
  [bonobo]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gtf.gz"
  [gibbon]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf.gz"
  [gorilla]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gtf.gz"
  [borangutan]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gtf.gz"
  [sorangutan]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gtf.gz"
)
for sp in "${!urls[@]}"; do
    mkdir -p "$sp"
    (
      cd "$sp"
      wget -q "${urls[$sp]}"
    ) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################################################################################
#Extract coding exons annotated in the GTF

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo $species
gunzip genes/"$species"/*.gtf.gz
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
#Extract coding exon coordinates from GTF
python3 scripts/gtf_to_bed.py genes/"$species"/"$gtf".gtf genes/"$species"/"$gtf".bed
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#########################################################################################################################
#bin the SCFR counts by length and GC content

cd /media/aswin/SCFR/SCFR-main
cd /media/aswin/SCFR/SCFR-main/SCFR_all
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
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
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Plot a 2 dimensional histogram of length vs AT content and label the SCFR longer than 10 Kb that overlap coding exons
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo $species
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
bedtools intersect -a SCFR_all/"$species"_long_SCFRs.bed -b genes/"$species"/"$gtf".bed -wa -wb|cut -f 5,13,17|sort -u > SCFR_all/"$species"_genes_of_interest.txt
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

cp scripts/combined_2dhist.r SCFR_all
cd SCFR_all
Rscript combined_2dhist.r

#########################################################################################################################
#Identify gene deserts 

mkdir /media/aswin/SCFR/SCFR-main/gene_deserts
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo $species
bed=$(find genes/$species/ -name "GCF_*.bed")
python3 scripts/gene_desert_finder.py $bed --z 2 --out gene_deserts/$species"_gene_deserts"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e










####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
####################################################################################################################################################################################################################################################################################################################

#Download large genome files
datasets download genome accession $genome --include genome,gtf,seq-report --dehydrated --filename $genome.zip
unzip $genome.zip -d $genome
time datasets rehydrate --directory $genome






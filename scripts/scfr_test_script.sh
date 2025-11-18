####################################################################################################################################################################################################################################################################################################################
#
####################################################################################################################################################################################################################################################################################################################

####################################################################################################################################################################################################################################################################################################################
#1. Prepare and collect data & script

#Download github repository
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
  mkdir -p QC/genome_metadata
  mkdir -p SCFR_lists
  gene_deserts
done

#mkdir github
#mv SCFR_lists PCA github/

#List of genomic accesions for analysis
echo -e "GCA_009914755.4\thuman\nGCA_028858775.3\tchimpanzee\nGCA_029289425.3\tbonobo\nGCA_028878055.3\tgibbon\nGCA_029281585.3\tgorilla\nGCA_028885625.3\tborangutan\nGCA_028885655.3\tsorangutan" > QC/genome_accessions

####################################################################################################################################################################################################################################################################################################################
#2. Download all the seven primate genomes

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
#3. Quality check

cd /media/aswin/SCFR/SCFR-main
mkdir -p QC/genome_metadata

#Get genome metadata using datasets
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
#Get length of chromosomes 
cd /media/aswin/SCFR/SCFR-main/QC/genome_metadata
find . -name "*_genome_metadata.txt" | xargs -n1 sh -c 'awk "{print\$9,\$8}" $0' | grep -v "refseq_accession" | sed '1i refseq_accession Length' | grep -v "^\-" > all_species_chr_length

#Calculate length of chromosome from fetched genome sequences
cd /media/aswin/SCFR/SCFR-main/chrs
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
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

#Report errors in input: Compare lengths between fetched and database chromosomes
cd /media/aswin/SCFR/SCFR-main/QC/genome_metadata
for chr in $(awk '{print$1}' all_species_chr_length | egrep -v "refseq_accession|-")
     do
     m1=$(awk -v a="$chr" '$1==a' all_species_chr_length)
     m2=$(awk -v a="$chr" '$2==a' fetched_length)
     echo $m1 $m2
     unset m1 m2
done | awk '{if($2==$5) print$0,"same"; else print$0,"different"}' |  sed '1i ncbi_accession ncbi_length species fetched_accession fetched_length Compare_length' | column -t > compare_chromosome_length

####################################################################################################################################################################################################################################################################################################################
#4.Identify and summarize SCFRs

#4.1. Identify SCFRs (51.6833 mins)
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

#4.2. Collate all gaps
cd /media/aswin/SCFR/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  echo $species
  cat SCFR/gaps/"$species"/gap*.bed|sort -k1,1 -k2n,2 > SCFR_all/gaps_"$species".bed
done

#4.3.1. Get summary (35.4167 mins)
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

#4.3.2. Since running 6 species need huge memory run 3 at a time (1.67528 hours ~100.517 mins)
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
#5. Strand assymetry

#5.1. Get strand assymetry of SCFRs per chromosome (11.3 mins)
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

#5.2. Get strand assymetry of SCFRs in sliding windows across the genomes (9.45 mins)
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

#5.3. Get summary & plot (6 secs)
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
#6. get GC content of the SCFRs (48.1667)

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
#7. Download all the seven primate genome annotation files

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
#8. Extract coding exons annotated in the GTF

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
#9. Bin the SCFR counts by length and GC content

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
#10. Gene deserts

#10.1. Identify gene deserts (1 min)
mkdir /media/aswin/SCFR/SCFR-main/gene_deserts
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo $species
bed=$(find genes/$species/ -name "GCF_*.bed")
python3 scripts/gene_desert_finder.py $bed --z 2 --out gene_deserts/$species"_intronic_intergenic"
python3 my_scripts/gene_desert_finder_intergenic.py $bed --z 2 --out gene_deserts/$species"_only_intergenic"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#------------------------------------------------------------------------------------------------------------------------------------------------------
#10.2. Plot length stats of gene deserts

#Get genome sizes
mkdir -p genome_sizes
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
    out=genome_sizes/${species}.genome
    > "$out"
    for fa in chrs/$species/*.fasta
    do
        chr=$(basename "$fa" .fasta)
        len=$(awk '/^>/ {next} {L+=length($0)} END {print L}' "$fa")
        echo -e "${chr}\t${len}" >> "$out"
    done
    echo "Created $out"
done

#Get gene desert bed file
for f in gene_deserts/*_only_intergenic_gene_deserts.tsv; do
    species=$(basename "$f" | cut -d'_' -f1)
    awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3}' "$f" > gene_deserts/${species}_only_intergenic_gene_deserts.bed
done

#Get SCFR in bed file
for f in SCFR_all/*SCFR_all.out; do
    species=$(basename "$f" | cut -d'_' -f1)
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$f" > SCFR_all/${species}_SCFR_all.bed
done

#Plot data
cd /media/aswin/SCFR/SCFR-main/
python3 compute_desert_stats.py gene_deserts/*.bed
Rscript plot_gene_desert_stats_2.r all_desert_lengths.tsv

#------------------------------------------------------------------------------------------------------------------------------------------------------
#10.3. Fishers test

#Merge SCFRs into genome-wide intervals
time sort -k1,1 -k2,2n SCFR_all/human_SCFR_all.out > human_SCFR_all_sorted.out
time bedtools merge -i human_SCFR_all_sorted.out | awk '{print$1,$2,$3}' OFS="\t" > human_SCFR_all_sorted_merged.bed
#Keep only long SCFRs
time awk '$3-$2 >= 1000' human_SCFR_all_sorted_merged.bed > human_SCFR_all_sorted_merged_atleast_1kb.bed
#Fisher's test
bedtools fisher -a human_SCFR_all_sorted_merged_atleast_1kb.bed -b gene_deserts/human_only_intergenic_gene_deserts.bed -g genome_sizes/human.genome


####################################################################################################################################################################################################################################################################################################################
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




####################################################################################################################################################################################################################################################################################################################
#12. Quantification of codon usage patterns and PCA:

#Create SCFR beds with different cut-off lengths

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
for win in 5000 7500 10000
do
echo " -"$win
#For SCFRs of atleast a specific length
mkdir -p PCA/"$species"/"$win"/with_coding_region
for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | cut -f 1|sort -u`
do
cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ > PCA/"$species"/"$win"/with_coding_region/"$chr".fasta
done
#For SCFRs of atleast a specific length and with no overlap with coding region
mkdir -p PCA/"$species"/"$win"/without_coding_region
for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | cut -f 1|sort -u`
do
cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ > PCA/"$species"/"$win"/without_coding_region/"$chr".fasta
done
#Calculate codon metrics
python3 scripts/codon_usage_metrics.py PCA/${species}/${win}/with_coding_region PCA/${species}/${win}/with_coding_region
python3 scripts/codon_usage_metrics.py PCA/${species}/${win}/without_coding_region PCA/${species}/${win}/without_coding_region
Rscript plotPCA.r PCA/"$species"/$win/with_coding_region PCA/"$species"/$win/with_coding_region
Rscript plotPCA.r PCA/"$species"/$win/without_coding_region PCA/"$species"/$win/without_coding_region
done
done
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






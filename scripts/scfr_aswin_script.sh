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
  mkdir -p Fourier_analysis
  mkdir -p gene_deserts
  mkdir -p SCFR_summaries
done

#mkdir github
#mv SCFR_lists PCA mv SCFR_summaries github/

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

#get genome metadata
cd /media/aswin/SCFR/SCFR-main/QC
awk '{print$1}' genome_accessions \
| xargs -n1 sh -c 'datasets summary genome accession $0 --as-json-lines | json2xml | xtract -pattern root -def "-" -element organism_name assembly_name assembly_level accession contig_l50 contig_n50 scaffold_l50 scaffold_n50 total_number_of_chromosomes total_sequence_length total_ungapped_length genome_coverage' \
| sed 's/ /_/g' | sed '1i Species Assembly_name Assembly_level Genbank_accession Refseq_accession Contig_l50 Contig_n50 Scaffold_l50 Scaffold_n50 Total_number_of_chromosomes Total_sequence_length Total_ungapped_length Genome_coverage' \
| sed 's/PRJNA[^ \t]\+//g' | sed 's/SAMN[^ \t]\+//g' | awk '!($1=="Homo_sapiens" || $1=="Species") { $1=""; sub(/^ +/, ""); }1' > genome_metadata/genome_metadata

#Get chromosome metadata using datasets
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

#Combine the summaries
ls | egrep "_SCFR_all_frame_summary_by_chromosome|_SCFR_all_frame_summary_all" | xargs -n1 sh -c 'mv $0 SCFR_summaries/'
cd /media/aswin/SCFR/SCFR-main/SCFR_summaries
ls | grep "_SCFR_all_frame_summary_all.csv" | xargs -n1 sh -c 'grep -v Chromosome $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_all_frame_summary_all.csv//g' | tr "," " " | tr -d '"' | sed '1i Species Chromosome Frame N Min Q1 Median Mean Q3 Max SD P95 P99 Q_1Kb Q_5Kb Q_10Kb' > all_species_chromsome_frame_wise_summary
ls | grep "_SCFR_all_frame_summary_all.csv" | xargs -n1 sh -c 'grep ALL $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_all_frame_summary_all.csv//g' | tr "," " " | tr -d '"' | sed '1i Species Chromosome Frame N Min Q1 Median Mean Q3 Max SD P95 P99 Q_1Kb Q_5Kb Q_10Kb' > all_species_frame_wise_summary

#Make SCFR length stats plot
#NOTE: Rename species common names before plotting
cd /media/aswin/SCFR/SCFR-main/SCFR_all
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/plot_all_species_SCFR_length_stats.R all_scfr_length_stats.tsv all_scfr_length_stats.pdf

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

#combine the summaries
cd /media/aswin/SCFR/SCFR-main/SCFR_all
ls | grep "_SCFR_asymmetries_out.csv" | xargs -n1 sh -c 'grep -v chrom $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_asymmetries_out.csv//g' | tr "," " " | tr -d '"' \
	| sed '1i Species chrom f_count r_count strand_count_asym f_length r_length strand_length_asym frame1_count frame1_length frame2_count frame2_length frame3_count frame3_length frame-1_count frame-1_length frame-2_count frame-2_length frame-3_count frame-3_length' > all_species_chromosome_wise_SCFR_asymmetries
ls | grep "_SCFR_asymmetries_out.csv" | xargs -n1 sh -c 'grep ALL $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_asymmetries_out.csv//g' | tr "," " " | tr -d '"' \
	| sed '1i Species chrom f_count r_count strand_count_asym f_length r_length strand_length_asym frame1_count frame1_length frame2_count frame2_length frame3_count frame3_length frame-1_count frame-1_length frame-2_count frame-2_length frame-3_count frame-3_length' > all_species_SCFR_asymmetries

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
#SCFRs from both intronic & intergenic region
python3 scripts/gene_desert_finder.py $bed --z 2 --out gene_deserts/$species"_intronic_intergenic"
#SCFRs from only intergenic region
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

cd /media/aswin/SCFR/SCFR-main
/media/aswin/SCFR/SCFR-main/gene_deserts/fishers_test
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  (
  echo " -"$species
  mkdir /media/aswin/SCFR/SCFR-main/gene_deserts/fishers_test/$species
  for tsv in $(find gene_deserts/ -maxdepth 1 -mindepth 1 -name "$species*.tsv")
  do
  prefix=$(echo $tsv | sed 's/.tsv//g' | awk -F "/" '{print$NF}')
  awk '!/start/ {print$1,$2,$3,$4}' OFS="\t" $tsv > gene_deserts/fishers_test/$species/$prefix".bed"
  /media/aswin/programs/bedtools2-2.31.1/bin/bedtools fisher -a human_SCFR_atleast_10000.out -b human_only_intergenic_intergenic_regions.bed -g /media/aswin/SCFR/SCFR-main/genome_sizes/human.genome
  

####################################################################################################################################################################################################################################################################################################################
#11. Filter SCFRs that are part of coding region

#Filter with different length cut offs (16.0167 mins)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for win in 0 100 500 1000 2500 5000 7500 10000
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
#Get window-wise summary of SCFR count and length & it's overlap with whole genome & cds from all species (99.3667 mins)

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
genome_size=$(awk '{a+=$2} END{print a}' /media/aswin/SCFR/SCFR-main/genome_sizes/$species".genome")
#cds
cd /media/aswin/SCFR/SCFR-main/genes/$species
cds=$(ls | grep "GC.*.bed")
#merge cds
	sort -k1,1V -k2,2n $cds | /media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i - > $species"_cds_merged.bed"
merged_cds_len=$(awk '{print$3-$2}' $species"_cds_merged.bed" | awk '{a+=$1} END{print a}')
#scfr
total_unfiltered_scfr_len=$(wc -l /media/aswin/SCFR/SCFR-main/SCFR_lists/0/$species"_SCFR_atleast_0.out" | awk '{print$1}')
cd /media/aswin/SCFR/SCFR-main/SCFR_lists
for win in 0 100 500 1000 2500 5000 7500 10000
do
cd $win
#get stats
	all_scfr=$(ls | grep $species"_SCFR_atleast_"$win".out")
	noncoding_scfr=$(ls | grep $species"_SCFR_atleast_"$win"_in_non_coding.bed")
	all_count=$(wc -l $all_scfr | awk '{print$1}')
	noncoding_count=$(wc -l $noncoding_scfr | awk '{print$1}')
	all_len=$(awk '{print$3-$2}' $all_scfr | awk '{a+=$1} END{print a}')
	noncoding_len=$(awk '{print$3-$2}' $noncoding_scfr | awk '{a+=$1} END{print a}')
	all_merged_len=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i $all_scfr | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
	noncoding_merged_len=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i $noncoding_scfr | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
#merge csfr
	sort -k1,1V -k2,2n /media/aswin/SCFR/SCFR-main/SCFR_lists/"$win"/"$all_scfr" | /media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i - > /media/aswin/SCFR/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed"
#get percentages
	per_unfiltered_scfr_by_filtered=$(calc $all_count / $total_unfiltered_scfr_len | tr -d "~" | awk '{print$1*100}')
	per_genome_by_scfr=$(calc $all_merged_len / $genome_size | tr -d "~" | awk '{print$1*100}')
	merged_scfr_len=$(awk '{print$3-$2}' /media/aswin/SCFR/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed" | awk '{a+=$1} END{print a}')
	scfr_cds_overlap=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a /media/aswin/SCFR/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed" -b /media/aswin/SCFR/SCFR-main/genes/$species/$species"_cds_merged.bed" | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
	per_scfr_by_coding=$(calc $scfr_cds_overlap / $all_merged_len | tr -d "~" | awk '{print$1*100}')
scfr=$(echo $species $win $all_count $noncoding_count $all_len $noncoding_len $all_merged_len $noncoding_merged_len $per_unfiltered_scfr_by_filtered $genome_size $per_genome_by_scfr $merged_cds_len $scfr_cds_overlap $per_scfr_by_coding)
echo $scfr
unset all_scfr noncoding_scfr all_count noncoding_count all_len noncoding_len all_merged_len noncoding_merged_len per_unfiltered_scfr_by_filtered per_genome_by_scfr merged_scfr_len scfr_cds_overlap per_scfr_by_coding 
cd /media/aswin/SCFR/SCFR-main/SCFR_lists
done
unset genome_size win cds merged_cds_len total_unfiltered_scfr_len
cd /media/aswin/SCFR/SCFR-main
done | sed '1i Species Window Total_No_SCFR Total_No_noncoding_SCFR Total_length_SCFR Total_length_noncoding_SCFR Merged_length_SCFR Merged_length_noncoding_SCFR Percent_unfiltered_by_filtered_SCFR Genome_size Percent_genome_by_SCFR Merged_CDS_length SCFR_CDS_overlap Percent_SCFR_by_coding' > /media/aswin/SCFR/SCFR-main/SCFR_summaries/window_wise_all_species_scfr_coding_stats_with_zero_window
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e > /media/aswin/SCFR/SCFR-main/SCFR_summaries/runtime_window_wise_all_species_scfr_coding_stats

#Plot SCFR count & length stats & percent in genome & cds 
awk -F "," '{print$1,$2,(($3-$4)/$3*100),($4/$3*100),$3,$9,$11,$14}' window_wise_all_species_scfr_coding_stats.csv | awk 'BEGIN{FS=OFS=","} NR==1{$3="Percent_coding_SCFR"; $4="Percent_noncoding_SCFR"}1' | tr " " "," > filtered_window_wise_all_species_scfr_coding_stats.csv
awk -F "," 'NR>1{print$1,$2,$3,(($3-$4)/$3*100),($4/$3*100),$9,$11,$14}' window_wise_all_species_scfr_coding_stats.csv | sed '1i Species,Window,Total_No_SCFR,Percent_coding_SCFR_count,Percent_noncoding_SCFR_count,Percent_unfiltered_by_filtered_SCFR,Percent_genome_by_SCFR,Percent_SCFR_by_coding' | sed 's/ /,/g' > filtered_window_wise_all_species_scfr_coding_stats.csv

#Get coding gene length stats of species
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
gene_length_stats=$(myfasta -l $species/GC*_cds.fa | awk '{print$NF}' | ministat -n | tail -1 | sed 's/^x //g' | sed 's/[ ]\+/ /g')
echo $species $gene_length_stats
unset species gene_length_stats
done | column -t > gene_length_stats


####################################################################################################################################################################################################################################################################################################################
#12. Quantification of codon usage patterns and PCA (10.7167 mins)

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
for win in 5000 7500 10000
do
echo " -"$win
#For SCFRs of atleast a specific length
mkdir -p PCA/"$species"/"$win"/with_coding_region
for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | cut -f 1|sort -u`
do
cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,1,"-"; else print$0,1,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > PCA/"$species"/"$win"/with_coding_region/"$chr".fasta
done
#For SCFRs of atleast a specific length and with no overlap with coding region
mkdir -p PCA/"$species"/"$win"/without_coding_region
for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | cut -f 1|sort -u`
do
cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,"-"; else print$0,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > PCA/"$species"/"$win"/without_coding_region/"$chr".fasta
done
#Calculate codon metrics
python3 scripts/codon_usage_metrics.py PCA/${species}/${win}/with_coding_region PCA/${species}/${win}/with_coding_region
python3 scripts/codon_usage_metrics.py PCA/${species}/${win}/without_coding_region PCA/${species}/${win}/without_coding_region
Rscript scripts/plotPCA.r PCA/"$species"/$win/with_coding_region PCA/"$species"/$win/with_coding_region
Rscript scripts/plotPCA.r PCA/"$species"/$win/without_coding_region PCA/"$species"/$win/without_coding_region
done
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#-----------------------------------------------------------------------------------------------------------
#For easier visual inspection plot PCA without labels (2.46667 mins)

mkdir /media/aswin/SCFR/SCFR-main/PCA_without_labels
cd /media/aswin/SCFR/SCFR-main/
cp -r PCA/* PCA_without_labels/
find PCA_without_labels/ -name "*.pdf" | xargs rm
find PCA_without_labels/ -name "*.fasta" | xargs rm

start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	for win in 5000 7500 10000
	do
	Rscript my_scripts/plotPCA_without_labels.r PCA_without_labels/"$species"/$win/with_coding_region PCA_without_labels/"$species"/$win/with_coding_region
	Rscript my_scripts/plotPCA_without_labels.r PCA_without_labels/"$species"/$win/without_coding_region PCA_without_labels/"$species"/$win/without_coding_region
	done
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

####################################################################################################################################################################################################################################################################################################################
#13. Discrete Fourier Transform Analysis

#-----------------------------------------------------------------------------------------------------------
#For SCFRs (for 5000, 7500 & 10000 - 72.8667 mins) (for 1000 & 2500 - )

cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
for win in 1000 2500 5000 7500 10000
do
echo " -"$win
#For SCFRs of atleast a specific length
	mkdir -p Fourier_analysis/"$species"/"$win"/with_coding_region
	for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | cut -f 1|sort -u`
	do
	cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,1,"-"; else print$0,1,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > Fourier_analysis/"$species"/"$win"/with_coding_region/"$chr".fasta
	done
#For SCFRs of atleast a specific length and with no overlap with coding region
	mkdir -p Fourier_analysis/"$species"/"$win"/without_coding_region
	for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | cut -f 1|sort -u`
	do
	cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,"-"; else print$0,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > Fourier_analysis/"$species"/"$win"/without_coding_region/"$chr".fasta
	done
#Run Fourier analysis
	for type in with_coding_region without_coding_region
	do
	cd Fourier_analysis/"$species"/"$win"/"$type"/
	for scfr in $(find . -mindepth 1 -maxdepth 1 -name "*.fasta" | cut -f2 -d "/")
	do
	echo " -"$scfr
	python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o "output_"$scfr -t 32 $scfr
	done
	unset scfr
	cd /media/aswin/SCFR/SCFR-main
	done
unset chr type
done
unset win
cd /media/aswin/SCFR/SCFR-main
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#-----------------------------------------------------------------------------------------------------------
#For genes

#Get GTF associated with genome
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes
cat /media/aswin/SCFR/SCFR-main/QC/genome_accessions | xargs -n2 bash -c 'paste <(echo $1 $0) <(datasets summary genome accession $0 --as-json-lines | json2xml | xtract -pattern root -def "-" -element accession)' | column -t > genome_gtf_accessions
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species

#Download cds from ncbi (21.9667 mins)
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes
start_time=$(date +%s)
while read i
	do
	species=$(echo $i | awk '{print$1}')
	gtf=$(echo $i | awk '{print$NF}')
	echo ">"$species $gtf
	time datasets download genome accession $gtf --include cds --filename $gtf".zip"
	unzip -o $gtf".zip" -d cds
	mkdir $species
	mv cds/ncbi_dataset/data/$gtf/cds_from_genomic.fna $species/$gtf"_cds.fa"
	unset species gtf
done < <(grep borangutan genome_gtf_accessions)
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

rm -r cds

#Run fourier analysis of genes (474.25)
cd /media/aswin/SCFR/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	cd Fourier_analysis/genes/$species
	gene=$(ls | grep "GC.*_cds.fa")
	python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o "output_"$gene -t 32 $gene
	cd /media/aswin/SCFR/SCFR-main
	unset gene
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#human genes failed to finish & no pdf was generated, likely due to large number of fasta in single file
#Split human cds fasta into 2 equal parts & run fourier seperately on both
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/split_human
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/split_human
cp ../human/GCF_009914755.1_cds.fa .
seqkit split -p 2 GCF_009914755.1_cds.fa
cd GCF_009914755.1_cds.fa.split/
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o output_GCF_009914755.1_cds.part_001.fa -t 32 GCF_009914755.1_cds.part_001.fa









####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
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



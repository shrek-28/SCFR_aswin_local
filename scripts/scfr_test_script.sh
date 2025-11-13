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
  time datasets download genome accession $genome --include genome,gtf,seq-report --filename $genome".zip"
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
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
     do
     cd $species
     for chr in $(ls *.fasta)
     do
     length=$(myfasta -l $chr | awk '{print$1,$NF}')
     echo $species $ch $length
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
#DRAFT SCRIPTS

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

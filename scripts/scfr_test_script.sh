


cd /media/aswin/SCFR/SCFR-main

time while read species
do
genome=$(echo $species | awk '{print$1}')
name=$(echo $species | awk '{print$2}')
echo ">"$genome $name
#Download
time datasets download genome accession $genome --include genome,gtf,seq-report --filename $genome".zip"
unset genome name
done < QC/genome_accessions


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



#
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

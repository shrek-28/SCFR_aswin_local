

#Download
time datasets download genome accession GCA_028878055.3 --include genome,gtf,seq-report --filename GCA_028878055.3.zip

unzip GCA_028878055.3.zip
cd ncbi_dataset/data/GCA_028878055.3

#Get accessions
cat sequence_report.jsonl | json2xml | xtract -pattern opt -def "-"  -element genbankAccession refseqAccession > genbank_refseq_accession

#Add refseq accesion in genome
awk 'NR==FNR {map[$1]=$2; next}
     /^>/ {
         hdr=$0
         for (k in map) {
             if (hdr ~ ">"k"($|[[:space:]])") {
                 sub(">"k, ">"map[k]" "k)
                 break
             }
         }
     }
     {print}' genbank_refseq_accession GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna \
> updated.fna



#

time while read chr 
do
g=$(echo $chr | awk '{print$1}')
r=$(echo $chr | awk '{print$2}')
echo $r $g
sed "/$g/ s/^>/>$r /g" GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna -i
done < genbank_refseq_accession

#

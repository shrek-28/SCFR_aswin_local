
cd /media/aswin/SCFR/SCFR-main/SCFR_summaries

mkdir /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot


cd /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot
find /media/aswin/SCFR/SCFR-main/genes/$species/ -name "*merged*" | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/SCFR_summaries/circos_plot/'
report=$(readlink -f /media/aswin/SCFR/SCFR-main/genome_reports/* | grep $species)
awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1'  $report | sed 's/ /_/g' | awk 'NR>1{print$9,$12}' | sed '1i accession names' > $species"_chromosome_names"
cp /media/aswin/SCFR/SCFR-main/genome_sizes/$species".genome" .


##########


find /media/aswin/SCFR/SCFR-main/SCFR_lists/ -name "human*" | egrep -v "in_non_coding" | xargs -n1 sh -c 'cp $0 .'
cp /media/aswin/SCFR/SCFR-main/genes/human/human_cds_merged.bed


find /media/aswin/SCFR/SCFR-main/SCFR_lists/ -name "human*" | egrep -v "in_non_coding" 
/media/aswin/SCFR/SCFR-main/genome_sizes/human.genome
/media/aswin/SCFR/SCFR-main/genes/human/human_cds_merged.bed



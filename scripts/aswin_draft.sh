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
awk -F "," -v a="$i" '$2==a {print$1,$3}' window_wise_all_species_scfr_coding_stats.csv | grep $s | awk '{print$2}'
done)
echo $i $j
unset j
done | sed '1i Window human bonobo chimpanzee gorilla borangutan sorangutan gibbon' | awk 'NR==1 {print; next} {for(i=2;i<=NF;i++) $i=$i/1000000; print}' OFS="\t" | column -t

#View higher strand asymmetry
sed 's/ /_/g' genome_reports/GCA_029289425.3_bonobo.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '{print$9,$3,$4,$5,$11}' | sort -k5,5n | colnum.sh 



cd /media/aswin/SCFR/SCFR-main/SCFR_summaries
find /media/aswin/SCFR/SCFR-main/SCFR_lists/ -name "human*" | egrep -v "in_non_coding|_SCFR_atleast_0.out" | xargs -n1 sh -c 'cp $0 circos_plot/'
cp /media/aswin/SCFR/SCFR-main/genome_sizes/human.genome circos_plot/
cp /media/aswin/SCFR/SCFR-main/genes/human/human_cds_merged.bed circos_plot/



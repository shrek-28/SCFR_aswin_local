##Plotting presence of Gene desert and intergenic regions in chhoromosome of one species (eg bonobo)
#Load datasets 
chrom_lengths <- read.table("C:/Users/dipan/IISER_project/SCFR_clone-main/genome_sizes/bonobo.genome",
                            +                             header=FALSE, sep="\t")
colnames(chrom_lengths) <- c("chr","length")
 chrom_lengths <- chrom_lengths %>% mutate(length_Mb = length/1e6)
 
gd <- read.table(file.choose(), header=TRUE, sep="\t") %>% mutate(type="Gene Desert")

ir <- read.table(file.choose(), header=TRUE, sep="\t") %>% mutate(type="Intergenic Region")
# Combine
> segments <- bind_rows(gd, ir) %>%
  +     mutate(start_Mb = start/1e6, end_Mb = end/1e6)
> # --- Box positions ---
  > types <- c("Gene Desert", "Intergenic Region")
> box_height <- 2       # thick boxes
> space <- 0.5
> y_positions <- data.frame()
> current_y <- 0
> 
  > for(chr in rev(chrom_lengths$chr)){
         for(t in types){ 
           y_positions <- rbind(y_positions, data.frame(chr=chr, type=t, ymin=current_y, ymax=current_y+box_height))
             current_y <- current_y + box_height + space
       }
    }

#Main CODE for bonobo
species_name <- "Bonobo"
ggplot() +
  geom_rect(data=chr_boxes, aes(xmin=0, xmax=length_Mb, ymin=ymin, ymax=ymax),
            fill="grey90", color="black", size=0.5) +
  geom_rect(data=bonobo_segments, aes(xmin=start_Mb, xmax=end_Mb, ymin=ymin, ymax=ymax, fill=type),
            color=NA) +
  scale_fill_manual(values=c("Gene Desert"="#FF6666", "Intergenic Region"="#66CCFF")) +
  scale_y_continuous(
    breaks = sapply(seq_along(unique(chr_boxes$chr)), function(i) {
      rows <- y_positions %>% filter(chr==unique(chr_boxes$chr)[i])
      mean(c(min(rows$ymin), max(rows$ymax)))
    }),
    labels = unique(chr_boxes$chr),
    expand = c(0,0)
  ) +
  scale_x_continuous(expand=c(0,0), limits=c(0, max(chr_boxes$length_Mb)*1.05)) +
  labs(x="Position (Mb)", y="Chromosome", fill="Region Type",
       title=paste0("Gene Desert & Intergenic Regions in ", species_name)) +
  theme_minimal() +
  theme(
    panel.grid=element_blank(),
    axis.text.y=element_text(face="bold"),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(face="bold"),
    legend.position="top",
    plot.title=element_text(face="bold", hjust=0.5)
  )

## Plotting presence of Gene desert and intergenic regions in chhoromosome for all seven primates

library(dplyr)
library(ggplot2)

species_list <- c("bonobo", "chimpanzee", "gorilla", "human", "borangutan", "gibbon", "sorangutan")

for(sp in species_list){
  
  # --- Genome lengths ---
  genome_file <- paste0("C:/Users/dipan/IISER_project/SCFR_clone-main/genome_sizes/", sp, ".genome")
  chrom_lengths <- read.table(genome_file, header=FALSE, sep="\t")
  colnames(chrom_lengths) <- c("chr","length")
  chrom_lengths <- chrom_lengths %>% mutate(length_Mb = length/1e6)
  
  # --- GD and IR files ---
  gd_file <- paste0("C:/Users/dipan/IISER_project/SCFR_clone-main/gene deserts/", sp, "_only_intergenic_gene_deserts.tsv")
  ir_file <- paste0("C:/Users/dipan/IISER_project/SCFR_clone-main/gene deserts/", sp, "_only_intergenic_intergenic_regions.tsv")
  
  gd <- read.csv(gd_file, sep="\t", header=TRUE) %>% mutate(type="Gene Desert")
  ir <- read.csv(ir_file, sep="\t", header=TRUE) %>% mutate(type="Intergenic Region")
  
  segments <- bind_rows(gd, ir) %>% mutate(start_Mb = start/1e6, end_Mb = end/1e6)
  
  # --- Box positions ---
  types <- c("Gene Desert", "Intergenic Region")
  box_height <- 2
  space <- 0.5
  y_positions <- data.frame()
  current_y <- 0
  
  for(chr in rev(chrom_lengths$chr)){
    for(t in types){
      y_positions <- rbind(y_positions, data.frame(chr=chr, type=t, ymin=current_y, ymax=current_y+box_height))
      current_y <- current_y + box_height + space
    }
  }
  
  # --- Filter NW chromosomes ---
  chrom_lengths$chr <- trimws(chrom_lengths$chr)
  segments$chr <- trimws(segments$chr)
  y_positions$chr <- trimws(y_positions$chr)
  
  chrom_lengths <- chrom_lengths %>% filter(!grepl("^NW", chr))
  segments <- segments %>% filter(!grepl("^NW", chr))
  y_positions <- y_positions %>% filter(!grepl("^NW", chr))
  
  # --- Merge y positions ---
  segments <- segments %>% left_join(y_positions, by=c("chr","type"))
  chr_boxes <- expand.grid(chr=chrom_lengths$chr, type=types) %>%
    left_join(y_positions, by=c("chr","type")) %>%
    left_join(chrom_lengths, by="chr") %>%
    mutate(xmin=0, xmax=length_Mb)
  
  # --- Plot ---
  p <- ggplot() +
    geom_rect(data=chr_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="grey90", color="black", size=0.5) +
    geom_rect(data=segments, aes(xmin=start_Mb, xmax=end_Mb, ymin=ymin, ymax=ymax, fill=type),
              color=NA) +
    scale_fill_manual(values=c("Gene Desert"="#FF6666","Intergenic Region"="#66CCFF")) +
    scale_y_continuous(
      breaks = sapply(seq_along(unique(chr_boxes$chr)), function(i) {
        rows <- y_positions %>% filter(chr==unique(chr_boxes$chr)[i])
        mean(c(min(rows$ymin), max(rows$ymax)))
      }),
      labels = unique(chr_boxes$chr),
      expand = c(0,0)
    ) +
    scale_x_continuous(expand=c(0,0), limits=c(0, max(chr_boxes$length_Mb)*1.05)) +
    labs(x="Position (Mb)", y="Chromosome", fill="Region Type",
         title=paste0("Gene Desert & Intergenic Regions in ", sp)) +
    theme_minimal() +
    theme(
      panel.grid=element_blank(),
      axis.text.y=element_text(face="bold"),
      axis.ticks.y=element_blank(),
      axis.text.x=element_text(face="bold"),
      legend.position="top",
      plot.title=element_text(face="bold", hjust=0.5)
    )
  
  # --- Save plot ---
  ggsave(filename=paste0("C:/Users/dipan/IISER_project/SCFR_clone-main/plots/", sp, "_gene_deserts.png"),
         plot=p, width=16, height=8)
  
  print(paste(sp, "plot done!"))
}




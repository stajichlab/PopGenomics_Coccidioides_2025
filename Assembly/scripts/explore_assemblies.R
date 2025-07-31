#!/usr/bin/env R
library(ggplot2)
library(tidyverse)
library(purrr)
library(readr)
library(cowplot) # it allows you to save figures in .png file
library(RColorBrewer)
#library(smplot2)

asminfo <- read_tsv("asm_stats.tsv",col_names=TRUE) %>% mutate(ID=str_replace(SampleID, ".AAFTF", ""))
samples <- read_csv("samples.csv",col_names=TRUE) 
merged <- asminfo %>% left_join(samples,by = join_by(ID == Strain))
unique(merged$Organism)
colourCount = length(unique(merged$Organism))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)

pdf("plots/size_explore_asm.pdf")
p <- ggplot(data = merged, mapping = aes(x = Organism, y = TOTAL_LENGTH/10**6,color=Organism)) +
  geom_boxplot() + 
  geom_jitter(width = 0.15) +  theme_cowplot(12) + theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_manual(values=getPalette) +
  ylab("Length (Mb)") 
p

ggsave("genome_size_boxplot.pdf",p,width=5,height=3)

p2 <- ggplot(data = merged %>% filter(BUSCO_Complete > 75), 
             mapping = aes(x = log(Average_Coverage)/log(10), y=BUSCO_Complete,color=Organism)) +
  geom_point()  +  theme_cowplot(12) + 
  scale_color_manual(values=getPalette) + xlab("Coverage (log(10))") + ylab("BUSCO Completeness") 

p2
ggsave("genome_coverage.pdf",p2,width=5,height=3)

p3 <- ggplot(data = merged %>% filter(BUSCO_Complete > 75), 
             mapping = aes(x = TOTAL_LENGTH/10**6, y=BUSCO_Duplicate,color=Organism)) +
  geom_point()  +  theme_cowplot(12)  + 
  scale_color_manual(values=getPalette) + xlab("Genome size (Mb)") + ylab("BUSCO Duplication")

p3
ggsave("genome_dup_by_size.pdf",p3,width=3,height=3)



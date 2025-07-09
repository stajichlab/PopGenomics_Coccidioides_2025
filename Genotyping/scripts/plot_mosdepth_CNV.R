#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(stringr)
library(cowplot)
library(purrr)
colors1 <- colorRampPalette(brewer.pal(8, "RdYlBu"))
manualColors = c("dodgerblue2", "red1", "grey20")

utils::globalVariables(c("pos", "Depth", "alternating_colors"))

alternating_colors = rep( c("red", "black", "blue", 
"orange","darkgreen","brown","slategray","purple"), times= 20)
Prefix = "GG70491" # this is Cimmitis prefix
mosdepthdir = "coverage/mosdepth"
ifelse( !dir.exists("results"),dir.create("results"),FALSE)
targetresultdir = file.path("results","CNV_plots")
ifelse( !dir.exists(targetresultdir),dir.create(targetresultdir),FALSE)


sampleinfo = read_csv('samples.csv',col_names=TRUE) %>% select(c(Strain,Organism)) %>% mutate(orgstr = sub(' ','_',Organism))

plot_strain <- function(strain, data) {
  l = subset(data, data$Strain == strain)
	Title = sprintf("Chr coverage plot for %s", strain)
	p <- ggplot(l, aes(x = pos, y = NormalizedDepth, color = CHR)) +
		scale_colour_manual(values = alternating_colors) +
		geom_point(alpha = 0.9, size = 0.8, shape = 16,show.legend = FALSE) +
		labs(title = Title,xlab = "Position", y = "Normalized Read Depth") +
		scale_x_continuous(name = "Chromosome",	expand = c(0, 0), breaks = ticks, labels = unique(l$CHR) ) +
		scale_y_continuous(name = "Normalized Read Depth",expand = c(0, 0), limits = c(0, 3)) +
		theme_classic()
  return(p)
		#+ guides(fill = guide_legend(keywidth = 3,keyheight = 1))
}

plot_chrs <- function(chrom, data) {
	Title = sprintf("Chr%s depth of coverage", chrom)
	l <- subset(data, data$CHR == chrom)
	l$bp <- l$Start
	p <- ggplot(l, aes(x = bp, y = NormalizedDepth, color = Strain)) +
	geom_point(alpha = 0.7,	size = 0.8, shape = 16) + # scale_color_brewer(palette='RdYlBu',type='seq') +
	labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
	scale_x_continuous(expand = c(0,	0), name = "Position") +
	scale_y_continuous(name = "Normalized Read Depth", expand = c(0, 0), limits = c(0, 3)) +
	theme_classic()
	#guides(fill = guide_legend(keywidth = 3,keyheight = 1))
}

orgs = unique(sampleinfo$orgstr)
#species = orgs[1]
for (species in orgs ) {
#for (window in windows) {
  window = 10000
  inpattern = sprintf(".%sbp.regions.bed.gz$", window)
  file_list_all <- list.files(path = mosdepthdir, pattern = inpattern)
  file_list = sampleinfo %>% filter(orgstr == species) %>% 
    mutate( fname=sprintf("%s.%sbp.regions.bed.gz",Strain,window)) %>% select(c(Strain,fname))
  
  file_list = file_list[file_list$fname %in% file_list_all,]
  
  flist <- file_list %>% mutate(fname=sprintf("%s/%s",mosdepthdir,fname))
  
  df <- flist$fname %>% map_dfr( read_tsv, .id = "source", col_names=c("Chr", "Start", "End", "Depth")) %>% 
    mutate(Strain = flist[source,]$Strain) %>% select(-c(source))
  
  mediandepth <- df %>%  group_by(Strain) %>%
    summarise(median_depth = median(Depth, na.rm = TRUE))
  
  bedwindows <- df %>% left_join(mediandepth) %>% mutate(NormalizedDepth = Depth / median_depth)
  
  length(bedwindows$Strain)
  length(unique(bedwindows$Strain))
  #colnames(bedwindows) = c("Chr", "Start", "End", "Depth", "Strain")
  bedwindows$CHR <- sub(Prefix, "", bedwindows$Chr, perl = TRUE)
  unique(bedwindows$CHR)
  bedwindows <- subset(bedwindows, as.numeric(bedwindows$CHR) < 7 )
  bedwindows$CHR <- factor(as.numeric(bedwindows$CHR))
  unique(bedwindows$CHR)
  bedwindows <- bedwindows[ with(bedwindows, order(CHR)), ]
  covsum <- bedwindows %>%  group_by(Strain,CHR) %>% 
      summarise(median_depth = median(NormalizedDepth),
                mean_depth = mean(NormalizedDepth),
                chrlen=max(End)) %>% 
  	  filter(!is.na(mean_depth))
  	# filter short chroms?
  	covsum_aneu <- covsum %>% select(CHR,Strain,mean_depth,chrlen) %>%  
  	  filter(mean_depth > 1.5 | mean_depth < 0.5) %>% 
  	  filter(chrlen > 100000 ) %>% 
  	  dplyr::mutate_if(is.numeric, round, 2) 
  	
  	cnvpath = file.path(targetresultdir,species)
  	ifelse( !dir.exists(cnvpath),dir.create(cnvpath),FALSE)
  	
  	write_tsv(covsum_aneu,sprintf("%s/strain_aneu_mean_coverage.txt",cnvpath))	
  	
  	covsum_aneu <- covsum %>% select(CHR,Strain,median_depth,chrlen) %>% 
    filter(median_depth > 1.5 | median_depth < 0.5) %>% 
    filter(chrlen > 100000 ) %>% 
    dplyr::mutate_if(is.numeric, round, 2)
  	
  	write_tsv(covsum_aneu,sprintf("%s/strain_aneu_median_coverage.txt",cnvpath))
  	Strains_With_Anneuploid = unique(covsum_aneu$Strain)	
  	
  	covsum_wide <- covsum %>%  select(CHR,Strain,mean_depth,chrlen) %>% 
      pivot_wider(names_from = Strain, values_from = mean_depth, values_fill = 0) %>% 
      arrange(as.numeric(CHR)) %>% dplyr::mutate_if(is.numeric, round, 2)
  	write_tsv(covsum_wide,sprintf("%s/strain_mean_coverage.txt",cnvpath))
  
  	covsum_wide <- covsum %>% select(CHR,Strain,median_depth,chrlen) %>% 
  	  pivot_wider(names_from = Strain, values_from = median_depth, values_fill = 0) %>% 
  	  arrange(as.numeric(CHR)) %>% dplyr::mutate_if(is.numeric, round, 2)
    write_tsv(covsum_wide,sprintf("%s/strain_median_coverage.txt",cnvpath))
  
  
    d = bedwindows
    d <- d[order(as.numeric(d$CHR), d$Start), ]
    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start, d$CHR, length))
    d$pos = NA
  
    nchr = length(unique(bedwindows$CHR))
    lastbase = 0
    ticks = NULL
    minor = vector(, 8)
  
    for (i in 1:nchr) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$Start
      } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
        lastbase = lastbase + max(d[d$index == (i - 1), "Start"])
        minor[i] = lastbase
        d[d$index == i, "Start"] = d[d$index == i, "Start"] - min(d[d$index ==
          i, "Start"]) + 1
        d[d$index == i, "End"] = lastbase
        d[d$index == i, "pos"] = d[d$index == i, "Start"] + lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    # ticks
    minorB <- tapply(d$End, d$index, max, probs = 0.5)
    # minorB minor
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
  
    #pdf(pdffile, width = 16, height = 4)
    Title = "Depth of sequence coverage"
    # let's remove some poor behaving strains w low mapping - could try to calc this from data but easier for now to just drop these 2
   #d = subset(d,! (d$Strain %in% c("CCFEE_5069","CCFEE_5042","CCFEE_5056","DBVPG_6094","EXF_10533")))
   Strainlist = unique(d$Strain)
   total = length(Strainlist)
   chunks = split(Strainlist, ceiling(seq_along(Strainlist)/40))
   for(i in 1:length(chunks)) {
    dsubStrains = d[d$Strain %in% chunks[[i]], ]
    theseticks <- tapply(dsubStrains$pos, dsubStrains$index, quantile, probs = 0.5)
     p <- ggplot(dsubStrains, aes(x = pos, y = NormalizedDepth, color = CHR)) + 
         geom_vline(mapping = NULL,xintercept = minorB, alpha = 0.5, linewidth = 0.1, colour = "grey15") +
         geom_point(alpha = 0.8, size = 0.4, shape = 16) +
         facet_wrap(~Strain,ncol=1) + scale_colour_manual(values = alternating_colors) +
         labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
         scale_x_continuous(name = "Chromosome",  expand = c(0, 0), breaks = theseticks, labels = unique(dsubStrains$CHR)) +
         scale_y_continuous(name = "Normalized Read Depth",expand = c(0, 0), limits = c(0, 3)) +
         theme_classic() + theme(legend.position="none")
     	pdffile = sprintf("%s/genomewide_coverage_%dkb.StrainSet%d.pdf", cnvpath,window/1000,i)
      ggsave(p,file=pdffile,width=20,height=30)
   }
  
   plts <- lapply(unique(d$Strain), plot_strain, data = d)
  
  	strains = unique(d$Strain)
  	for(i in 1:length(unique(d$Strain) ) ) {
  		pdffile=sprintf("%s/StrainPlot_%dkb.%s.pdf", cnvpath,window/1000,strains[[i]])
    	ggsave(plot = plts[[i]], file = pdffile, width=20, height=8)
  	}
  	
  	anneu = d[d$Strain %in% Strains_With_Anneuploid,]
  	StrainsAnn = unique(anneu$Strain)
    if ( length(StrainsAnn) > 0 ) {
  		theseticks <- tapply(anneu$pos, anneu$index, quantile, probs = 0.5)
  	  p <- ggplot(anneu, aes(x = pos, y = NormalizedDepth, color = CHR)) + 
  	  geom_vline(mapping = NULL,xintercept = minorB, alpha = 0.5, linewidth = 0.1, colour = "grey15") +
  	  geom_point(alpha = 0.8, size = 0.4, shape = 16) +
  	  facet_wrap(~Strain,ncol=1) + scale_colour_manual(values = alternating_colors) +
  	  labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
  	  scale_x_continuous(name = "Chromosome",  expand = c(0, 0), breaks = theseticks, labels = unique(anneu$CHR)) +
  	  scale_y_continuous(name = "Normalized Read Depth",expand = c(0, 0), limits = c(0, 3)) +
  	  theme_classic() + theme(legend.position="none")
  	  ggsave(sprintf("%s/Anneuploid_Strains_Wrap.pdf",cnvpath),p,width=20,height=30)
  	
  	  annplts <- lapply(StrainsAnn, plot_strain, data = anneu)
  	  for(i in 1:length(StrainsAnn ) ) {
  	    pdffile=sprintf("%s/Anneu_StrainPlot_%dkb.%s.pdf", cnvpath,window/1000,StrainsAnn[i])
    	  ggsave(plot = annplts[[i]], file = pdffile, width=20, height=8)
    	}
  }
	
}
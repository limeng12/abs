setwd("/Users/mengli/Documents/projects/abs");
library(valr)
library(dplyr)
source("code/multiplot.r")
source("hepg2/code/get_hepg2_high_exp_trans.R",echo=TRUE)

####high express gene should###
genomefile <- valr_example('hg19.chrom.sizes.gz');

genome <- read_genome(genomefile);

###only use intron in high expressed gene(TPM>10)  to avoid bias
#intron<-read_bed("anno/intron_coor_gencode_high_hepg2.bed",n_fields = 6);

intron<-read_bed("anno/hg19_gencode_intron_from_ucsc.bed",n_fields = 6);
intron<-intron[!duplicated(intron[,c("chrom","start","end")]  ),];
intron$start<-intron$start+1


intron<-intron[intron$name %in% high_exp_trans_ids,];



#intron<-intron[(intron$end-intron$start>250) & (intron$end-intron$start<50000),];
intron<-intron[(intron$end-intron$start>500),];


meta_5_3<-function(y_file){
  print(y_file)
  y<-read_bedgraph(paste0("/Volumes/mengli/abs/hepg2/eclip_bedgraph_hepg2/",y_file) );
  
  s5ss <- intron %>%
    filter(strand == '+') %>%
    mutate(end = start + 1)
  
  
  s3ss <- intron %>%
    filter(strand == '+') %>%
    mutate(start = end - 1)
  
  
  region_size <- 200
  # 50 bp windows
  win_size <- 2
  
  # add slop to the TSS, break into windows and add a group
  x5 <- s5ss %>%
    bed_slop(genome, both = region_size) %>%
    bed_makewindows(win_size);
  
  x3 <- s3ss %>%
    bed_slop(genome, both = region_size) %>%
    bed_makewindows(win_size);
  
  
  res5 <- bed_map(x5, y, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    summarize(win_mean = mean(win_sum, na.rm = TRUE),
              win_sd = sd(win_sum, na.rm = TRUE) );
  
  res3 <- bed_map(x3, y, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    summarize(win_mean = mean(win_sum, na.rm = TRUE),
              win_sd = sd(win_sum, na.rm = TRUE) );
  
  
  
  
  library(ggplot2)
  
  x5_labels <- seq(-region_size, region_size, by = win_size * 5);
  x5_breaks <- seq(1, max(x5$.win_id), by = 5);
  sd5_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)
  
  #x5_labels[x5_labels<0]<-str_c("exon_",x5_labels[x5_labels<0])
  
  #x5_labels[x5_labels>0]<-str_c("intron_",x5_labels[x5_labels>0])
  
  p5<-ggplot(res5, aes(x = .win_id, y = win_mean)) +
    geom_point() + #geom_pointrange(sd5_limits) + 
    scale_x_continuous(labels = x5_labels, breaks = x5_breaks) + 
    xlab('Position (bp from junction)') + ylab('Signal') + 
    ggtitle(paste0(y_file," exon-intron")) +
    theme_classic();
  
  
  x3_labels <- seq(-region_size, region_size, by = win_size * 5);
  x3_breaks <- seq(1, max(x3$.win_id), by = 5);
  sd3_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)
  
  p3<-ggplot(res3, aes(x = .win_id, y = win_mean)) +
    geom_point() + #geom_pointrange(sd3_limits) + 
    scale_x_continuous(labels = x3_labels, breaks = x3_breaks) + 
    xlab('Position (bp from junction)') + ylab('Signal') + 
    ggtitle(paste0(y_file," intron-exon")) +
    theme_classic();
  
  
  #print(p5);
  multiplot(p5,p3,cols=1)
  
}


pdf("hepg2/result/jc_meta_hepg2.pdf", width=8, height = 4);

files_all<-list.files("/Volumes/mengli/abs/hepg2/eclip_bedgraph_hepg2/");

#files_all<-c("AQR_1plus.bedgraph","AQR_2plus.bedgraph");
for(i in files_all){
  meta_5_3(i);
  #break
}

dev.off();


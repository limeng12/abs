#  control 
#
setwd("/Users/mengli/Documents/projects/abs");
library(valr)
library(plyr)
library(dplyr)
library(bedr)
library(stringr)
source("code/multiplot.r");
#source("code/jc53_peak_meta_minus.R",echo = TRUE);


system("sh fly/code/run_sh/generate_bedgraph_from_cryptic_ss_fly.sh");

####high express gene should###
# genomefile <- valr_example('hg19.chrom.sizes.gz');

genome <- read_genome("fly/anno/chrom_sizes.txt");


###only use intron in high expressed transcript (TPM>10) and in protein coding genes to avoid bias
#intron<-read_bed("anno/intron_coor_gencode_high.bed",n_fields = 6);
intron<-read_bed("fly/anno/dm6_ensembl_intron_from_ucsc.bed",n_fields = 6);

intron<-intron[!duplicated(intron[,c("chrom","start","end")]  ),];

intron$chrom<-str_sub(intron$chrom,4)

intron<-intron[!str_detect(intron$chrom,"_"),];


intron$start<-intron$start+1;

#intron<-intron[-1* which(is.na(intron$start)),];


intron<-intron[(intron$end-intron$start>500),];

#intron_regon<-str_c(intron$chrom,":",intron$start,"-",intron$end);

print(paste0("Number of bp intron used: ", nrow(intron) ) );



meta_5_3<-function(rbp){
  #  rbp<-"3_ENCSR015TOB_x16"
  #if(!str_detect(y_file,"plus.bedgraph") | str_detect(y_file,"input") ){
  #  return(0);
  #}
  
  #y_file<-"AQR_1plus.bedgraph";
  # rbp<-"TIA1_1"
  
  #"ZNF800_2plus.bedgraph"         
  #"ZRANB2_1minus.bedgraph"
  
  y<-read_bedgraph(paste0("/Users/mengli/Documents/projects/abs/fly/data/star_target_only_jc_bedgraph/",
                          rbp,"_no_ctl.bed.bedgraph") );
  print(rbp);
  
  s5ss <- intron %>%
    filter(strand == '+') %>%
    mutate(end = start );
  
  
  s3ss <- intron %>%
    filter(strand == '+') %>%
    mutate(start = end);
  
  
  
  region_size <-500
  # 50 bp windows
  win_size <- 5
  
  # add slop to the TSS, break into windows and add a group
  x5 <- s5ss %>%
    bed_slop(genome, left = 75,right=region_size) %>%
    bed_makewindows(win_size);
  
  x3 <- s3ss %>%
    bed_slop(genome, left=region_size,right=75) %>%
    bed_makewindows(win_size);
  
  
  
  res5 <- bed_map(x5, y, win_sum = sum(value, na.rm = TRUE) ) %>% 
    group_by(.win_id) %>% 
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE),
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  res3 <- bed_map(x3, y, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE),
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  
  
  ##run_minus(rbp);
  #y_minus<-read_bedgraph(paste0("/Volumes/mengli/abs/k562_eclip_coverage/",rbp,".bedgraph") );
  y_minus<-y;
  print(rbp)
  
  s5ss <- intron %>%
    filter(strand == '-') %>%
    mutate(start = end );
  
  
  s3ss <- intron %>%
    filter(strand == '-') %>%
    mutate(end = start);
  
  
  #region_size <- 300
  # 50 bp windows
  #win_size <- 1
  
  # add slop to the TSS, break into windows and add a group
  x5 <- s5ss %>%
    bed_slop(genome, left = region_size,right=75) %>%
    bed_makewindows(win_size);
  
  x3 <- s3ss %>%
    bed_slop(genome, left=75,right=region_size) %>%
    bed_makewindows(win_size);

  
  
  res5_minus <- bed_map(x5, y_minus, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE),
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  res3_minus <- bed_map(x3, y_minus, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE),
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  
  
  res5$win_sum<-res5$win_sum+rev(res5_minus$win_sum);
  
  res3$win_sum<-res3$win_sum+rev(res3_minus$win_sum);
  
  
  library(ggplot2);
  
  x5_labels <- seq(-75, region_size, by = win_size * 5);
  x5_breaks <- seq(0, max(x5$.win_id), by = 5);
  sd5_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd);
  
  #x5_labels[x5_labels<0]<-str_c("exon_",x5_labels[x5_labels<0]);
  
  #x5_labels[x5_labels>0]<-str_c("intron_",x5_labels[x5_labels>0]);
  
  p5<-ggplot(res5, aes(x = .win_id, y = win_sum) ) +
    geom_point() + geom_line()+#geom_pointrange(sd5_limits) + 
    scale_x_continuous(labels = x5_labels, breaks = x5_breaks) + 
    xlab('Position (bp from junction)') + ylab('Signal') + 
    ggtitle(paste0(rbp," exon-intron 5' -> 3' ") ) +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1) );
  
  
  x3_labels <- seq(-region_size, 75, by = win_size * 5);
  x3_breaks <- seq(0, max(x3$.win_id), by = 5);
  sd3_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)
  
  p3<-ggplot(res3, aes(x = .win_id, y = win_sum)) +
    geom_point() + geom_line()+#geom_pointrange(sd3_limits) + 
    scale_x_continuous(labels = x3_labels, breaks = x3_breaks) + 
    xlab('Position (bp from 3ss junction)') + ylab('Signal') + 
    ggtitle(paste0(rbp," intron-exon 5' -> 3' ")) +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1));
  
  
  #print(p5);
  multiplot(p5,p3,cols=1);
  
  rm(s5ss,s3ss,x5,x3,res5,res3);
  
}

pdf("fly/result/jc_meta_cryptic_site.pdf", width=15, height = 8);
library(stringr);

files_all<-list.files("/Users/mengli/Documents/projects/abs/fly/data/star_target_only_jc_bedgraph");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
gene_ids<-sapply(str_split(files_all,"\\_no_ctl"),"[",1);
#gene_ids<-unique(gene_ids[!str_detect(files_all,"input")]);
#gene_ids<-readLines("samples/gene_id.txt");
#rbp<-"ATP5C1"
#gene_ids<-c(str_c("3_",gene_pol_ids),str_c("5_",gene_pol_ids) );
#gene_ids<-gene_ids[1:1];
#gene_ids<-gene_ids[1:5];
#gene_ids<-"AQR";

for(i in gene_ids){
  meta_5_3(i);
  #break
}

dev.off();


#"ZNF800_2plus.bedgraph"         
#"ZRANB2_1minus.bedgraph"


setwd("/Users/mengli/Documents/projects/abs");
library(valr)
library(plyr)
library(dplyr)
library(bedr)
library(stringr)
source("code/multiplot.r");
#source("code/meta_profile/get_k562_high_exp_trans.R")
#source("code/jc53_peak_meta_minus.R",echo = TRUE);

system("sh code/run_sh/generate_bedgraph_from_cryptic_ss.sh");

####high express gene should###
genomefile <- valr_example('hg19.chrom.sizes.gz');

genome <- read_genome(genomefile);

###only use intron in high expressed transcript (TPM>10) and in protein coding genes to avoid bias
#intron<-read_bed("anno/intron_coor_gencode_high.bed",n_fields = 6);
intron<-read_bed("anno/hg19_gencode_intron_from_ucsc.bed",n_fields = 6);
#intron<-intron[!duplicated(intron[,c("chrom","start","end")]),]


intron$name<-sapply(str_split(intron$name,"\\."),"[",1);

##only consider exclusively express transcripts
#intron<-intron[intron$name %in% high_exp_trans_ids,];

intron<-intron[!duplicated(intron[,c("chrom","start","end")]  ),];

intron<-intron[!str_detect(intron$chrom,"_"),];


intron<-intron[(intron$end-intron$start>500),];

intron_regon<-str_c(intron$chrom,":",intron$start,"-",intron$end);

print(paste0("Number of bp intron used: ", nrow(intron) ) );


bp_site<-read_bed("anno/FileS2_BS_site.bed",n_fields = 6);

bp_site$score<-as.numeric(bp_site$score);

bp_site<- bp_site %>% filter(score > 1 );

bp_site_region<-str_c(bp_site$chrom,":",bp_site$start,"-",bp_site$end);


in_high_intron<-bp_site_region %in.region% intron_regon;


intron$start<-intron$start+1;

bp_site$start<-bp_site$start+1


bp_site<-bp_site[in_high_intron,];


print(paste0("Number of introns used: ", nrow(intron) ) );

print(paste0("Number of bp used: ", nrow(bp_site) ) );

#intron<-intron[1:80000,];
#bp_site<-bp_site[1:5000,];


meta_5_3<-function(rbp){
  
  
  y<-read_bedgraph(paste0("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/",
                          rbp,"_no_ctl.bed.bedgraph") );
  print(rbp);
  
  s5ss <- intron %>%
    filter(strand == '+') %>%
    mutate(end = start );
  
  
  s3ss <- intron %>%
    filter(strand == '+') %>%
    mutate(start = end);
  
  
  bp<- bp_site %>%
    filter(strand == '+') %>%
    mutate(start = end );
  
  
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
  
  xbp<-bp %>%
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
  
  resbp<-bed_map(xbp, y, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE) ,
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  y_minus<-y;
  print(rbp)
  
  s5ss <- intron %>%
    filter(strand == '-') %>%
    mutate(start = end );
  
  
  s3ss <- intron %>%
    filter(strand == '-') %>%
    mutate(end = start);
  
  
  bp<- bp_site %>%
    filter(strand == '-') %>%
    mutate(start = end );
  
  
  
  # add slop to the TSS, break into windows and add a group
  x5 <- s5ss %>%
    bed_slop(genome, left = region_size+1,right=75-1) %>%
    bed_makewindows(win_size);
  
  x3 <- s3ss %>%
    bed_slop(genome, left=75+1,right=region_size-1) %>%
    bed_makewindows(win_size);
  
  xbp<-bp %>%
    bed_slop(genome, left=75+1,right=region_size-1) %>%
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
  
  resbp_minus<-bed_map(xbp, y_minus, win_sum = sum(value, na.rm = TRUE) ) %>%
    group_by(.win_id) %>%
    dplyr::summarize(win_mean = mean(win_sum, na.rm = TRUE),
              #win_sd = sd(win_sum, na.rm = TRUE),
              win_sum = sum(win_sum, na.rm = TRUE)
    );
  
  
  
  res5$win_sum<-res5$win_sum+rev(res5_minus$win_sum);
  
  res3$win_sum<-res3$win_sum+rev(res3_minus$win_sum);
  
  resbp$win_sum<-resbp$win_sum+rev(resbp_minus$win_sum);
  
  
  
  
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
  
  
  xbp_labels <- seq(-region_size, 75, by = win_size * 5);
  xbp_breaks <- seq(0, max(xbp$.win_id), by = 5);
  sdbp_limits <- aes(ymax = win_mean + win_sd, ymin = win_mean - win_sd)
  
  pbp<-ggplot(resbp, aes(x = .win_id, y = win_sum)) +
    geom_point() +geom_line()+ #geom_pointrange(sd3_limits) + 
    scale_x_continuous(labels = xbp_labels, breaks = xbp_breaks) + 
    xlab('Position (bp from branch point)') + ylab('Signal') + 
    ggtitle(paste0(rbp," 5' -> 3'")) + 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1));
  
  
  #theme_classic();
  #print(p5);
  multiplot(p5,p3,pbp,cols=1);
  
 # rm(s5ss,s3ss,bp,x5,x3,xbp,res5,res3,resbp,y,y_minus);
  rm(s5ss,s3ss,bp,x5,x3,xbp,res5,res3,resbp,p5,p3,pbp,res5_minus,res3_minus,resbp_minus,y,y_minus);
  
}


unlink("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/3_rbp_merge.bedgraph")

files_all_3<-list.files("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph",
                        pattern="^3_*",full.names=TRUE);

cryptic_merge_3<-read.table(files_all_3[1],sep = "\t",as.is = TRUE,header = FALSE);
colnames(cryptic_merge_3)<-c("chr","start","end","count");

for(f in files_all_3[-1]){
  if(str_detect(f,"CTL") || str_detect(f,"R749H_slow") ||
     str_detect(f,"RBM8A_hela") || str_detect(f,"WT_poII") || str_detect(f,"rbp_merge")){
    next;
  }
  print(f)
  cryptic_merge_3_one<-read.table(f,sep = "\t",as.is = TRUE,header = FALSE)
  colnames(cryptic_merge_3_one)<-c("chr","start","end","count");
  
  
  cryptic_merge_3<-rbind(cryptic_merge_3, cryptic_merge_3_one);
  
  #break
}

cryptic_merge_3_one_gr<-cryptic_merge_3 %>% dplyr::group_by(chr,start,end) %>% dplyr::summarise(count=sum(count));


write.table(cryptic_merge_3_one_gr,col.names = FALSE,quote = FALSE,sep="\t",row.names = FALSE,
            file="/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/3_rbp_merge.bedgraph")



unlink("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/5_rbp_merge.bedgraph")
files_all_5<-list.files("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph",
                        pattern="^5_*",full.names=TRUE);

cryptic_merge_5<-read.table(files_all_5[1],sep = "\t",as.is = TRUE,header = FALSE);

cryptic_merge_5<-read.table(files_all_3[1],sep = "\t",as.is = TRUE,header = FALSE);
colnames(cryptic_merge_5)<-c("chr","start","end","count");


for(f in files_all_5[-1]){
  if(str_detect(f,"CTL") || str_detect(f,"R749H_slow") ||
     str_detect(f,"RBM8A_hela") || str_detect(f,"WT_poII") || str_detect(f,"rbp_merge")){
    next;
  }
  print(f)
  cryptic_merge_5_one<-read.table(f,sep = "\t",as.is = TRUE,header = FALSE)
  colnames(cryptic_merge_5_one)<-c("chr","start","end","count");
  
  
  cryptic_merge_5<-rbind(cryptic_merge_5, cryptic_merge_5_one);
  
  #break
}

cryptic_merge_5_one_gr<-cryptic_merge_5 %>% dplyr::group_by(chr,start,end) %>% dplyr::summarise(count=sum(count));


write.table(cryptic_merge_5_one_gr,col.names = FALSE,quote = FALSE,sep="\t",row.names = FALSE,
            file="/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/5_rbp_merge.bedgraph")

# 
# 
# pdf("result/jc_meta_cryptic_site_merge.pdf", width=15, height = 8);
# 
# y1<-read_bedgraph(paste0("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/3_rbp_merge.bedgraph") );
# 
# meta_5_3("3_all_RBP_merge", y1);
# 
# 
# y2<-read_bedgraph(paste0("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph/5_rbp_merge.bedgraph") );
# 
# meta_5_3("5_all_RBP_merge", y2);
# 
# 
# dev.off();


#"ZNF800_2plus.bedgraph"         
#"ZRANB2_1minus.bedgraph"


pdf("result/jc_meta_cryptic_site_raw.pdf", width=15, height = 8);
library(stringr);
files_all<-list.files("/Users/mengli/Documents/projects/abs/data/star_target_only_jc_bedgraph");
gene_ids_in_file<-sapply(str_split(files_all,"\\_no"),"[",1);

gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );
gene_ids<-c(gene_ids_bren,gene_ids_sicr);
gene_ids<-gene_ids[!str_detect(gene_ids,"_inte")]

select_genes<-c( str_c("3_", gene_ids), 
                 str_c("5_",gene_ids ) )


gene_ids<-intersect(gene_ids_in_file,select_genes);

for(i in gene_ids){
  if(str_detect(i,"CTL") ){
    next;
  }
  
  meta_5_3(i);
  #break
}

dev.off();



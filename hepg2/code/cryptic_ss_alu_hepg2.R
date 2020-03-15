library(plyr);
library(stringr);
library(dplyr);
#library(BSgenome);
#library(BSgenome.Hsapiens.UCSC.hg19);
library(ggplot2);
library(GenomicRanges)
library(IRanges)
options("scipen"=100, "digits"=4)
#hg19_genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19");
setwd("/Users/mengli/Documents/projects/abs");

gene_ids<-unique(readLines("hepg2/samples/gene_id_hepg2.txt") ) 

ctl_sj<-read.table("hepg2/data/star_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);

repeat_anno<-read.table("anno/hg19_repeat2",header = TRUE,as.is = TRUE,sep = "\t");

repeat_anno_alu<-repeat_anno[repeat_anno$repFamily=="Alu",];

repeat_anno_alu[,"length"]<-repeat_anno_alu$genoEnd-repeat_anno_alu$genoStart;

repeat_anno_alu_filter<-repeat_anno_alu[(repeat_anno_alu$length>250)|(repeat_anno_alu$length<350),];

repeat_anno_alu_filter_gr<-with(repeat_anno_alu_filter,
                                GRanges(seqnames = genoName,IRanges(start =genoStart ,end = genoEnd),strand=strand )  ) ;



gene_ids<-gene_ids[!str_detect(gene_ids,"CTL_")];

g_53_sj<-matrix(nrow=0,ncol=3);

for(g in gene_ids){
  print(g);
  
  #t_sj<-read.table(paste0("data/star_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  t_sj_5ss<-read.table(paste0("hepg2/data/star_target_only_jc_5ss/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE);
  colnames(t_sj_5ss)<-c("chr","start","end","X5_n","strand","V6");
  t_sj_5ss[,"start"]<-t_sj_5ss[,"start"]+1;
  
  
  t_sj_3ss<-read.table(paste0("hepg2/data/star_target_only_jc_3ss/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE);
  colnames(t_sj_3ss)<-c("chr","start","end","X3_n","strand","V6");
  t_sj_3ss[,"start"]<-t_sj_3ss[,"start"]+1;
  
  
  t_sj_gr_5ss<-with(t_sj_5ss, GRanges(seqnames = chr,IRanges(start=start,end=end) ,strand="*"))
  t_sj_gr_3ss<-with(t_sj_3ss, GRanges(seqnames = chr,IRanges(start=start,end=end) ,strand="*"))
  
  ss5_in_intron<-countOverlaps(t_sj_gr_5ss,repeat_anno_alu_filter_gr,type="within");
  ss3_in_intron<-countOverlaps(t_sj_gr_3ss,repeat_anno_alu_filter_gr,type="within");
  
  
  ss5_in_alu<-sum(ss5_in_intron>0);
  ss3_in_alu<-sum(ss3_in_intron>0);
  
  #ss5_not_in_exon_count_per<-ss5_not_in_exon_count/nrow(t_sj_5ss)
  #ss3_not_in_exon_count_per<-ss3_not_in_exon_count/nrow(t_sj_3ss)
  
  g_53_sj<-rbind(g_53_sj,c(g, ss5_in_alu, ss3_in_alu) );
  
}


g_53_sj_data<-as.data.frame(g_53_sj,stringsAsFactors=FALSE);
colnames(g_53_sj_data)<-c("rbp","cryptic_5ss_alu","cryptic_3ss_alu");


write.table(g_53_sj_data,file="hepg2/result/g_53_alu.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



gene_read_fr<-read.table("hepg2/result/rbp_read_count_map.tsv",sep = "\t",as.is = TRUE,header = TRUE);

g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr,by=c("rbp"="g") );

g_53_sj_data_read[,"normalized cryptic_5ss_alu"]<-as.numeric(g_53_sj_data_read[,"cryptic_5ss_alu"])/
  as.numeric(g_53_sj_data_read$read_count)*(10^6);

g_53_sj_data_read[,"normalized cryptic_3ss_alu"]<-as.numeric(g_53_sj_data_read[,"cryptic_3ss_alu"])/
  as.numeric(g_53_sj_data_read$read_count)*(10^6);

write.table(g_53_sj_data_read,file="hepg2/result/g_53_alu.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
library(ggplot2);
options("scipen"=100, "digits"=4)

gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

gene_ids<-c(gene_ids_bren,gene_ids_sicr);
gene_ids<-gene_ids[!str_detect(gene_ids,"_inte")]


g_53_sj<-matrix(nrow=0,ncol=2);

for(g in gene_ids){
  print(g);
  if(str_detect(g,"_inte") || str_detect(g,"_CTL_")){
    next;
  }
  
  jc_fr<-read.table( paste0("data/star_target_only_jc/",g,"_no_ctl.bed"),header = FALSE,as.is = TRUE );
  
  g_53_sj<-rbind(g_53_sj,c(g,sum(jc_fr[,8])) )
  
}
  
g_53_sj_data<-as.data.frame(g_53_sj,stringsAsFactors=FALSE);
colnames(g_53_sj_data)<-c("rbp","total_read_count")

gene_read_fr_fr<-read.table("result/rbp_read_count_map.tsv",header = TRUE,as.is = TRUE,sep="\t")

g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr_fr,by=c("rbp"="g") );

summary(as.numeric(g_53_sj_data_read[,"total_read_count"])/as.numeric(g_53_sj_data_read[,"read_count"]));



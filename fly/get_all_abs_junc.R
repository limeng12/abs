setwd("/Users/mengli/Documents/projects/abs");
options(scipen=200)
library(plyr);
library(stringr);
library(dplyr);

#gene_ids<-unique(readLines("fly/samples/gene_id_fly.txt") );
gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt"),readLines("fly/samples/gene_id_fly_all.txt") ) ));

#gene_ids<-unique(readLines("fly/samples/gene_id_fly_old.txt") );

ctl_sj<-read.table("fly/data/star_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);

#colnames(ctl_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno")


anno_sj<-read.table("fly/anno/dm6_ensembl_intron_from_ucsc.bed",sep = "\t",header = TRUE, as.is = TRUE);
anno_sj[,2]<-anno_sj[,2]+1;

gene_ids_noctl<-gene_ids[!str_detect(gene_ids,"CTL")];


all_junc<-matrix(nrow = 0,ncol=5);
  
abs3ss<-0;

for(g in gene_ids_noctl){
  print(g);
  
  t_sj<-read.table(paste0("fly/data/star_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  
  #colnames(t_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno");
  t_sj<-t_sj[(t_sj$anno>3) & (t_sj$type>0),];
  
  
  t_ctl_sj<-left_join(t_sj,ctl_sj,
                      by=c("chr"="chr","X5_pos"="X5_pos","X3_pos"="X3_pos",
                           "strand"="strand","type"="type","X5_pos"="X5_pos","X3_pos"="X3_pos") );
  
  t_only_sj_all<-t_ctl_sj[is.na(t_ctl_sj$anno.y),];
  
  abs3ss<-c(abs3ss,str_c(t_only_sj_all$chr,":",t_only_sj_all$X3_pos) );
  
  #####reverse the junction, since left is 5', right is 3' for junc
  is_neg<-t_only_sj_all[,"strand"]==2;
  
  a<-t_only_sj_all[is_neg,"X5_pos"];
  t_only_sj_all[is_neg,"X5_pos"]<-t_only_sj_all[is_neg,"X3_pos"];
  t_only_sj_all[is_neg,"X3_pos"]<-a;
  
  
  all_junc<-rbind(all_junc,cbind(t_only_sj_all[,1:4], rep(g,nrow(t_only_sj_all)))  );
  
  
}

all_junc_data<-as.data.frame(all_junc,stringsAsFactors=FALSE);
colnames(all_junc_data)<-c("chr","X5_pos","X3_pos","strand","rbp");

all_junc_nodup<-ddply(all_junc_data,.(chr,X5_pos,X3_pos,strand),function(x){paste0(x[,"rbp"],collapse = ",")} );

colnames(all_junc_nodup)<-c("chr","start","end","strand","rbp");

rbp_num<-str_count(all_junc_nodup[,"rbp"],",")+1;

all_junc_nodup<-cbind(all_junc_nodup,rbp_num);

is_positive<-all_junc_nodup[,"strand"]==1;
is_negative<-all_junc_nodup[,"strand"]==2;

all_junc_nodup[is_positive,"strand"]<-"+";
all_junc_nodup[is_negative,"strand"]<-"-";

cat(paste0("# of cryptic 3ss: ",length(unique(abs3ss)),"\n" ), file="result/count_of_cryptic_3ss.tsv",sep="\t" );

write.table(all_junc_nodup,file="fly/result/all_junc_nodup.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);


all_junc_nodup_bed<-all_junc_nodup;

all_junc_nodup_bed$start<-all_junc_nodup_bed$start-1;

score<-rep(0,nrow(all_junc_nodup_bed) );

all_junc_nodup_bed<-cbind(all_junc_nodup_bed,score);


all_junc_nodup_bed<-all_junc_nodup_bed[,c("chr","start","end","rbp","score","strand")];

#all_junc_nodup_bed_s<-all_junc_nodup_bed[1128000:(1000000+130000),];

write.table(all_junc_nodup_bed,file="fly/result/all_junc_nodup.bed",sep="\t",
            quote = FALSE,col.names = FALSE,row.names = FALSE);



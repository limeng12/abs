library(plyr);
library(stringr);
library(dplyr);
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
library(ggplot2);
library(GenomicRanges)
library(IRanges)
options("scipen"=100, "digits"=4)
#hg19_genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19");
setwd("/Users/mengli/Documents/projects/abs");

#gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
#gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

#gene_ids<-c(gene_ids_bren,gene_ids_sicr);

gene_ids<-unique(readLines("hepg2/samples/gene_id_hepg2.txt") ) 

ctl_sj<-read.table("hepg2/data/star_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);
#hg19_refgene.gtf


# setdiff(x, y, ignore.strand=FALSE)

exon_anno<-read.table("anno/gencode.v29lift37.annotation_exon.gtf",header = FALSE,as.is = TRUE,sep="\t");
#exon_anno<-read.table("anno/gencode.v29lift37.annotation_exon.gtf",header = FALSE,as.is = TRUE,sep="\t");

exon_anno<-gtf_anno[gtf_anno[,3]=="exon",]
#enst_ids<-sapply(str_split(exon_anno[,9],";"),"[",2);
#exon_anno[,"enst_ids_cut"]<-str_sub(enst_ids,16,30);

exon_anno<-exon_anno[,c(1,4,5)];
colnames(exon_anno)<-c("chr","start","end");
#exon_anno_sole_exp_sel<-unique(exon_anno_sole_exp_sel);
exon_anno_gr<-with(exon_anno,GRanges(seqnames = chr,IRanges(start=start,end=end) ,strand="*"))

#deep_intron_anno<-setdiff(trans_anno_gr,exon_anno_gr);


gene_ids<-gene_ids[!str_detect(gene_ids,"CTL_")];

g_53_sj<-matrix(nrow=0,ncol=5);

for(g in gene_ids){
  print(g);
  
  #t_sj<-read.table(paste0("data/star_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  t_sj_5ss<-read.table(paste0("hepg2/data/star_target_only_jc_5ss/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE);
  colnames(t_sj_5ss)<-c("chr","start","end","X5_n","strand","V6");
  t_sj_5ss[,"start"]<-t_sj_5ss[,"start"]+1;
  
  
  t_sj_3ss<-read.table(paste0("hepg2/data/star_target_only_jc_3ss/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE);
  colnames(t_sj_3ss)<-c("chr","start","end","X3_n","strand","V6");
  t_sj_3ss[,"start"]<-t_sj_3ss[,"start"]+1;
  
  
  t_sj_5ss$strand<-sapply(t_sj_5ss$strand,function(x){
    if(x==1){
      return("+");
    }
    return("-");
  });
  t_sj_3ss$strand<-sapply(t_sj_3ss$strand,function(x){
    if(x==1){
      return("+");
    }
    return("-");
  });
  
  
  t_sj_gr_5ss<-with(t_sj_5ss, GRanges(seqnames = chr,IRanges(start=start,end=end) ,strand="*"))
  t_sj_gr_3ss<-with(t_sj_3ss, GRanges(seqnames = chr,IRanges(start=start,end=end) ,strand="*"))
  
  ss5_in_exon<-countOverlaps(t_sj_gr_5ss,exon_anno_gr,type="within");
  ss3_in_exon<-countOverlaps(t_sj_gr_3ss,exon_anno_gr,type="within");
  
  
  ss5_in_exon_count<-sum(ss5_in_exon>0);
  ss3_in_exon_count<-sum(ss3_in_exon>0);
  
  
  if(ss5_in_exon_count>0){
    
    t_5_seqs<-getSeq(hg19_genome,t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="+"),"chr"],
                     t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="+"),"start"]-3,
                     t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="+"),"start"]+5,
                     strand=t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="+"),"strand"],as.character=TRUE);
    #cat(as.character(t_5_seqs),file=paste0("hepg2/data/star_abs5_motif_exon/",g,".seq"), sep="\n");
    
    t_5_seqs_neg<-getSeq(hg19_genome,t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="-"),"chr"],
                         t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="-"),"start"]-5,
                         t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="-"),"start"]+3,
                         strand=t_sj_5ss[(ss5_in_exon>0)&(t_sj_5ss$strand=="-"),"strand"],as.character=TRUE);
    #cat(as.character(t_5_seqs_neg),file=paste0("hepg2/data/star_abs5_motif_exon/",g,".seq"), sep="\n",append = TRUE);
    t_5_seqs<-c(t_5_seqs,t_5_seqs_neg)
    
    weblogo(t_5_seqs,file.out=paste0("hepg2/data/star_abs5_motif_exon/",g,".pdf") ,format="pdf" ,open = FALSE )
    
  }
  
  
  if(ss3_in_exon_count>0){
    t_3_seqs<-getSeq(hg19_genome,t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="+"),"chr"],
                     t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="+"),"start"]-19,
                     t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="+"),"start"]+3,
                     strand=t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="+"),"strand"],as.character=TRUE);
    #cat(as.character(t_3_seqs),file=paste0("hepg2/data/star_abs3_motif_exon/",g,".seq"), sep="\n");
    
    t_3_seqs_neg<-getSeq(hg19_genome,t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="-"),"chr"],
                         t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="-"),"start"]-3,
                         t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="-"),"start"]+19,
                         strand=t_sj_3ss[(ss3_in_exon>0)&(t_sj_3ss$strand=="-"),"strand"],as.character=TRUE);
    #cat(as.character(t_3_seqs_neg),file=paste0("hepg2/data/star_abs3_motif_exon/",g,".seq"), sep="\n",append = TRUE);
    
    
    t_3_seqs<-c(t_3_seqs,t_3_seqs_neg)
    weblogo(t_3_seqs,file.out=paste0("hepg2/data/star_abs3_motif_exon/",g,".pdf") ,format="pdf" ,open = FALSE )
    
  }
  
  
  
  
  ss5_in_exon_count_per<-ss5_in_exon_count/nrow(t_sj_5ss)
  ss3_in_exon_count_per<-ss3_in_exon_count/nrow(t_sj_3ss)
  
  
  g_53_sj<-rbind(g_53_sj,c(g,ss5_in_exon_count,ss3_in_exon_count,
                           ss5_in_exon_count_per,ss3_in_exon_count_per) );
  
  
}


g_53_sj_data<-as.data.frame(g_53_sj,stringsAsFactors=FALSE);
colnames(g_53_sj_data)<-c("rbp","cryptic_5ss_exon","cryptic_3ss_exon",
                          "cryptic_5ss_exon/total 5ss","cryptic_3ss_exon/total 3ss");

#g_53_sj_data[,"cryptic_5ss_over_3ss_not_near_exon"]<-as.numeric(g_53_sj_data[,"cryptic_5ss_not_near_exon"])/
#  as.numeric(g_53_sj_data[,"cryptic_3ss_not_near_exon"])

write.table(g_53_sj_data,file="hepg2/result/g_53_in_exon_full.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



gene_read_fr<-read.table("hepg2/result/rbp_read_count_map.tsv",sep = "\t",as.is = TRUE,header = TRUE);

g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr,by=c("rbp"="g") );

g_53_sj_data_read[,"normalized cryptic_5ss_exon"]<-as.numeric(g_53_sj_data_read[,"cryptic_5ss_exon"])/
  as.numeric(g_53_sj_data_read$read_count)*(10^6);

g_53_sj_data_read[,"normalized cryptic_3ss_exon"]<-as.numeric(g_53_sj_data_read[,"cryptic_3ss_exon"])/
  as.numeric(g_53_sj_data_read$read_count)*(10^6);

write.table(g_53_sj_data_read,file="hepg2/result/g_53_in_exon_full.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



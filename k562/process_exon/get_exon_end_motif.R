setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19);
library(RWebLogo)
options(scipen=999)

library(BSgenome.Hsapiens.UCSC.hg19);
genome <- BSgenome.Hsapiens.UCSC.hg19;


intron_coor<-read.table(file="anno/hg19_gencode_intron_from_ucsc.bed",sep="\t",
                        as.is = TRUE,header = FALSE);

intron_coor[,2]<-intron_coor[,2]+1;

colnames(intron_coor)<-c("chr","start","end","name","score","strand");

intron_coor_pos<-intron_coor[intron_coor[,6]=="+",]
intron_coor_range_positive<-with(intron_coor_pos,
                        GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                                strand=strand) );

intron_coor_neg<-intron_coor[intron_coor[,6]=="-",]
intron_coor_range_negative<-with(intron_coor_neg,
                                 GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                                         strand=strand) );




gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

gene_ids<-c(gene_ids_bren,gene_ids_sicr);

#gene_ids<-setdiff(gene_ids,"CTL_ENCSR000AEO");
#gene_ids<-setdiff(gene_ids,"ENCSR000KYM_FLAG_RNA");
#gene_ids<-setdiff(gene_ids,"ENCSR443QFD_CTL_RNA_RNAi_inte");





for(g in gene_ids){
  print(g);
  
  if(!file.exists(paste0("data/exon_mer_target_only_small_intron/",g,"_no_ctl.bed"))){
    print(paste0("don't have ",g))
    next;
  }
  
  t_sj<-read.table(paste0("data/exon_mer_target_only_intron/",g,"_no_ctl.bed"),sep = "\t",header = TRUE, as.is = TRUE);
  colnames(t_sj)<-c("chr","start","end","read_count","x1")
  
  t_only_sj_all_bed_small_range<-with(t_sj,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                                                                      strand="*",read_count=read_count) );
  
  in_positive_intron<-countOverlaps(t_only_sj_all_bed_small_range,intron_coor_range_positive,type="within");
  
  t_sj_positive_intron<-t_sj[in_positive_intron & (t_sj[,3]-t_sj[,2] >10 ),];
  
  if(nrow(t_sj_positive_intron)==0){
    next;
  }
  
  intron_seq_all_pos<-as.character(getSeq(genome,t_sj_positive_intron[,1],
                                          t_sj_positive_intron[,3]-10,t_sj_positive_intron[,3]+5+5 ,strand="+"));
  
  
  in_negative_intron<-countOverlaps(t_only_sj_all_bed_small_range,intron_coor_range_negative,type="within");
  
  t_sj_negative_intron<-t_sj[in_negative_intron & (t_sj[,3]-t_sj[,2] >10 ),];
  
  if(nrow(t_sj_negative_intron)==0){
    next;
  }
  
  intron_seq_all_neg<-as.character(getSeq(genome,t_sj_negative_intron[,1],
                                          t_sj_negative_intron[,2]-4-5,t_sj_negative_intron[,2]+11,strand="-"));
  
  intron_seq_all<-c(intron_seq_all_pos,intron_seq_all_neg)
  
  #intron_seq_all_fa<-str_c(">",1:length(intron_seq_all),"\n",intron_seq_all,"\n");
  
  #cat(intron_seq_all_fa,sep="",file=paste0("data/exon_mer_target_only_small_intron_fa/",g,".fa"))
   
  
  weblogo(intron_seq_all,file.out=paste0("data/exon_mer_target_only_small_intron_fa/",g,".pdf") ,open = FALSE )
  
}
  



#pdf("result/test_motif.pdf")

#dev.off();







  
  
  
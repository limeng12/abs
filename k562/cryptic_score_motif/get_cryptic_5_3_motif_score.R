## Weblogo motif of all cryptic site for each shRNA-seq

setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);
#library(BSgenome)
#library(BSgenome.Hsapiens.UCSC.hg19);
library(RWebLogo)
options(scipen=999);

#library(BSgenome.Hsapiens.UCSC.hg19);
#genome <- BSgenome.Hsapiens.UCSC.hg19;


#gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

#gene_ids<-c(gene_ids_bren,gene_ids_sicr);
gene_ids<-gene_ids_bren

#gene_ids<-setdiff(gene_ids,"CTL_ENCSR000AEO");
#gene_ids<-setdiff(gene_ids,"ENCSR000KYM_FLAG_RNA");
#gene_ids<-setdiff(gene_ids,"ENCSR443QFD_CTL_RNA_RNAi_inte");
#

#    gene_ids<-"ENCSR624OUI_AQR_poly"

g_53_sj<-matrix(nrow=0,ncol=3);

for(g in gene_ids){
  print(g);
  
  
  if(!file.exists(paste0("data/star_abs5_motif/",g,".seq"))){
    print(paste0("don't have ",g))
    next;
  }
  
  intron_seq_all<-readLines(paste0("data/star_abs5_motif/",g,".seq") );
  
  
  weblogo(intron_seq_all,file.out=paste0("data/star_abs5_motif/",g,".png") ,format="png" ,open = FALSE )
  weblogo(intron_seq_all,file.out=paste0("data/star_abs5_motif/",g,".pdf") ,format="pdf" ,open = FALSE )
  
  
  setwd("/Users/mengli/Documents/software/maxentscan")
  system(paste0("perl score5.pl /Users/mengli/Documents/projects/abs/data/star_abs5_motif/",g,".seq > /Users/mengli/Documents/projects/abs/data/star_abs5_motif/",g,".seq.score") )
  setwd("/Users/mengli/Documents/projects/abs");
  
  
  value_fr<-read.table(paste0("/Users/mengli/Documents/projects/abs/data/star_abs5_motif/",g,".seq.score"),sep="\t",header = TRUE,as.is = TRUE)[,2];
  
  median_score_5ss<-median(value_fr);
  
  
  
  if(!file.exists(paste0("data/star_abs3_motif/",g,".seq"))){
    print(paste0("don't have ",g))
    next;
  }
  
  intron_seq_all<-readLines(paste0("data/star_abs3_motif/",g,".seq"));
  #colnames(t_sj)<-c("chr","start","end","read_count","x1")
  
  weblogo(intron_seq_all,file.out=paste0("data/star_abs3_motif/",g,".png"),format="png" ,open = FALSE )
  weblogo(intron_seq_all,file.out=paste0("data/star_abs3_motif/",g,".pdf"),format="pdf" ,open = FALSE )
  
  
  setwd("/Users/mengli/Documents/software/maxentscan")
  system(paste0("perl score3.pl /Users/mengli/Documents/projects/abs/data/star_abs3_motif/",g,".seq > /Users/mengli/Documents/projects/abs/data/star_abs3_motif/",g,".seq.score") )
  setwd("/Users/mengli/Documents/projects/abs");
  
  value_fr<-read.table(paste0("/Users/mengli/Documents/projects/abs/data/star_abs3_motif/",g,".seq.score"),sep="\t",header = TRUE,as.is = TRUE)[,2];
  
  median_score_3ss<-median(value_fr);
  
  
  g_53_sj<-rbind(g_53_sj,c(g,median_score_5ss,median_score_3ss) );
  
  
}


g_53_sj_data<-data.frame(rbp=g_53_sj[,1],max_ent_5ss=g_53_sj[,2],max_ent_3ss=g_53_sj[,3] );
#library("plotflow")
write.table(g_53_sj_data,file="result/g_53_sj_data_max_ent_score.tsv", 
            quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE);




setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);

#gene_ids<-unique(readLines("fly/samples/gene_id_fly.txt") );
#gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt"),readLines("fly/samples/gene_id_fly_all.txt") ) ));
gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt") ) ));


files_all<-list.files("fly/data/star_log_fly/");

gene_read_fr<-matrix(nrow=0,ncol=2);

for( g in gene_ids){
  print(g);
  
  files_one_g<-files_all[str_detect(files_all,g)];
  
  star_log<-read.table(paste0("fly/data/star_log_fly/",files_one_g[1]),
                             sep = "\t",header = FALSE,fill = TRUE,as.is = TRUE);
  
  num_mapped_reads<-as.numeric(star_log[8,2]);
  
  for(i_file in files_one_g[-1]){
    
    star_log_t2<-read.table(paste0("fly/data/star_log_fly/",i_file),
                            sep = "\t",header = FALSE,fill = TRUE,as.is = TRUE);
    
    num_mapped_reads<-num_mapped_reads+as.numeric(star_log_t2[8,2]);
  }
  
  gene_read_fr<-rbind(gene_read_fr,c(g,num_mapped_reads));
  
}

gene_read_fr_fr<-as.data.frame(gene_read_fr,stringsAsFactors=FALSE)
colnames(gene_read_fr_fr)<-c("g","read_count");

write.table(gene_read_fr_fr,file="fly/result/rbp_read_count_map.tsv",
             sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)



  
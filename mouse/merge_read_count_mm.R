setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);


#source("code/process_star/get_star_read_count.R",echo=TRUE);

#   /Users/mengli/Documents/projects/abs/hepg2/data/star_log
#   find star_fly/ -type f -name '*.final.out' -exec cp '{}' star_fly_log/ ';'
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/fly/star_fly_log/\* /Users/mengli/Documents/projects/abs/fly/data/star_log_fly/

gene_ids<-unique(readLines("mouse/gene_id_mm.txt") );


files_all<-list.files("mouse/star_log/");

gene_read_fr<-matrix(nrow=0,ncol=3);

for( g in gene_ids){
  print(g);

  
  files_one_g<-files_all[str_detect(files_all,fixed(g))];
  
  if(length(files_one_g)==0){
    next;
  }
  
  star_log<-read.table(paste0("mouse/star_log/",files_one_g[1]),
                       sep = "\t",header = FALSE,fill = TRUE,as.is = TRUE);
  
  ave_read_length<-as.numeric(star_log[6,2]);
  #print(ave_read_length);
  num_mapped_reads<-as.numeric(star_log[8,2]);
  
  for(i_file in files_one_g[-1]){
    
    star_log_t2<-read.table(paste0("mouse/star_log/",i_file),
                            sep = "\t",header = FALSE,fill = TRUE,as.is = TRUE);
    
    num_mapped_reads<-num_mapped_reads+as.numeric(star_log_t2[8,2]);
  }
  
  gene_read_fr<-rbind(gene_read_fr,c(g,num_mapped_reads,ave_read_length));
  
}
            

gene_read_fr_fr<-as.data.frame(gene_read_fr,stringsAsFactors=FALSE);
colnames(gene_read_fr_fr)<-c("g","read_count","read_length");


#gene_read_fr_fr_raw<-read.table("mouse/rbp_read_count_map.tsv",header = TRUE,as.is = TRUE,sep="\t");

#gene_read_fr_fr<-rbind(gene_read_fr_fr,gene_read_fr_fr_raw);

write.table(gene_read_fr_fr,file="mouse/rbp_read_count_map.tsv",
            sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)




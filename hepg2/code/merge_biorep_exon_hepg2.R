setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);
library(bedr);
options(scipen=999)

#options("scipen"=100, "digits"=4);

#system("scp limeng@10.10.118.191:/picb/rnasys2/limeng/data/exon/\\*exon /Users/mengli/Documents/projects/abs/data/exon/")

region_len<-function(x){
  region_list<-str_split(x,":|-");
  start<-as.numeric(sapply(region_list,"[",2));
  end<-as.numeric(sapply(region_list,"[",3));
  len<-end-start+1;
  return(len);
}

gene_ids<-(unique(readLines("hepg2/samples/gene_id_hepg2.txt") ) );

flag_id<-"asdfaf"

# size.region

files_all<-list.files("hepg2/data/exon/");

for( g in gene_ids){
  print(g);
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  #sj_all<-read.table(paste0("data/exon/",files_one_g[1]),header = FALSE,as.is = TRUE)
  sj_all<-data.frame(chr=c(),start=c(),end=c(),id=c(),score=c(),strand=c() );

  for(i_file in files_one_g){
    if(file.size(paste0("hepg2/data/exon/",i_file))==0 ){
      next;
    }
    
    files_all_one_gene<-read.table(paste0("hepg2/data/exon/",i_file),header = FALSE,as.is = TRUE)
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
  colnames(sj_all)<-c("region","read_count");
  
  sj_all_mer<-sj_all %>% dplyr::group_by(region) %>% dplyr::summarise(read_count=sum(read_count))

  
  write.table(sj_all_mer,file = paste0("hepg2/data/exon_mer/",g,"_tab"),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = TRUE );
    
}


################################################merge all ctl#################################################

files_all<-list.files("hepg2/data/exon/");
sj_all<-data.frame(chr=c(),start=c(),end=c(),id=c(),score=c(),strand=c())

for( g in gene_ids){
  print(g);
  
  if(!str_detect(g,"CTL_") ){
    next;
  }
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  
  for(i_file in files_one_g){
    
    files_all_one_gene<-read.table(paste0("hepg2/data/exon/",i_file),header = FALSE,as.is = TRUE)
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
}  

colnames(sj_all)<-c("region","read_count");
  
sj_all_mer<-sj_all %>% dplyr::group_by(region) %>% dplyr::summarise(read_count=sum(read_count) );

write.table(sj_all_mer,file = paste0("hepg2/data/exon_mer/CTL_all_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );



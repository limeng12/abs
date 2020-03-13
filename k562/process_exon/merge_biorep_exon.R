setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);
library(bedr);
options(scipen=999)

#options("scipen"=100, "digits"=4);

#  scp limeng@10.10.118.191:/picb/rnasys2/limeng/data/exon/\*exon /Users/mengli/Documents/projects/abs/data/exon/


#  x<-sj_all_mer$region

region_len<-function(x){
  region_list<-str_split(x,":|-");
  start<-as.numeric(sapply(region_list,"[",2));
  end<-as.numeric(sapply(region_list,"[",3));
  len<-end-start+1;
  return(len);
}


#gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
#gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

#gene_ids<-c(gene_ids_bren,gene_ids_sicr);

#gene_ids<-setdiff(gene_ids,"CTL_ENCSR000AEO");
#gene_ids<-setdiff(gene_ids,"ENCSR000KYM_FLAG_RNA");
#gene_ids<-setdiff(gene_ids,"ENCSR443QFD_CTL_RNA_RNAi_inte");


flag_id<-"asdfaf"

# size.region

files_all<-list.files("hek293/exon/");

for( g in gene_ids){
  print(g);
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  #sj_all<-read.table(paste0("data/exon/",files_one_g[1]),header = FALSE,as.is = TRUE)
  sj_all<-data.frame(chr=c(),start=c(),end=c(),id=c(),score=c(),strand=c() );
  
  # sj_all<-read.table(paste0("data/circular/",files_one_g[1]),header = FALSE,as.is = TRUE)
  
  
  for(i_file in files_one_g){
    if(file.size(paste0("hek293/exon/",i_file))==0 ){
      next;
    }
    
    files_all_one_gene<-read.table(paste0("hek293/exon/",i_file),header = FALSE,as.is = TRUE)
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
  colnames(sj_all)<-c("region","read_count");
  
  
  # sj_1<-read.table(paste0("data/exon/",g,"_1_exon"),
  #                  header = FALSE,as.is = TRUE);
  # colnames(sj_1)<-c("region","read_count");
  # 
  # 
  # sj_2<-read.table(paste0("data/exon/",g,"_2_exon"),
  #                  header = FALSE,as.is = TRUE);
  # colnames(sj_2)<-c("region","read_count");
  #  
  #  
  #  
  # sj_all<-rbind(sj_1,sj_2);
  
  
  #sj_all_mer<-plyr::count(sj_all, vars = c("region") );
  sj_all_mer<-sj_all %>% dplyr::group_by(region) %>% dplyr::summarise(read_count=sum(read_count))
  
  #sj_all_mer<-ddply(sj_all, .(region),function(x){sum(x[,"read_count"]) } );
  #colnames(sj_all_mer)[2]<-"read_count" 
  
  #sj_len<-region_len(sj_all_mer[,"region"]);
  #sj_all_mer[sj_len<50,];
  
  write.table(sj_all_mer,file = paste0("hek293/exon_mer/",g,"_tab"),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = TRUE );
  
  #rm(sj_1,sj_1_filter,sj_2,sj_2_filter,sj_all,sj_all_mer);
  
}


################################################merge all ctl#################################################

# gene_ids_ctl<-gene_ids[str_detect(gene_ids,"CTL_")];
# 
# #gene_ids_ctl<-c(str_c(gene_ids_ctl,"_1") ,str_c(gene_ids_ctl,"_2") );
# #gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_norm_4");
# #""
#

files_all<-list.files("hek293/exon/");
sj_all<-data.frame(chr=c(),start=c(),end=c(),id=c(),score=c(),strand=c())

for( g in gene_ids){
  print(g);
  
  if(!str_detect(g,"CTL_") ){
    next;
  }
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  #sj_all<-read.table(paste0("data/exon/",files_one_g[1]),header = FALSE,as.is = TRUE)
  
  # sj_all<-read.table(paste0("data/exon/",files_one_g[1]),header = FALSE,as.is = TRUE)
  
  
  for(i_file in files_one_g){
    
    files_all_one_gene<-read.table(paste0("hek293/exon/",i_file),header = FALSE,as.is = TRUE)
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
}  


colnames(sj_all)<-c("region","read_count");
  
#sj_all_mer<-ddply(sj_all, .(region),function(x){sum(x[,"read_count"]) } );
  
#sj_all_mer<-plyr::count(sj_all, vars = c("region") );
#colnames(sj_all_mer)[2]<-"read_count" 

sj_all_mer<-sj_all %>% dplyr::group_by(region) %>% dplyr::summarise(read_count=sum(read_count) );


#sj_len<-region_len(sj_all_mer[,"region"]);

#sj_all_mer[sj_len<50,];

#colnames(sj_all_mer)[2]<-"read_count";


write.table(sj_all_mer,file = paste0("hek293/exon_mer/CTL_all_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );



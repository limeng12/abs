setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);

#   find star/ -type f -name '*SJ.out.tab' -exec cp '{}' star_re/ ';'
#   find star/ -type f -name '*.final.out' -exec cp '{}' star_log/ ';'

#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/hepg2/star_re/\*SJ.out.tab /Users/mengli/Documents/projects/abs/hepg2/data/star/
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/hepg2/star_log/\*.final.out /Users/mengli/Documents/projects/abs/hepg2/data/star_log/

gene_ids<-unique(readLines("hepg2/samples/gene_id_hepg2.txt") ); 

flag_id<-"nothing"


files_all<-list.files("hepg2/data/star/");

for( g in gene_ids){
  print(g);
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  if(length(files_one_g)==0){
    print(str_c("fail: ",g));
    next;
  }
  
  sj_all<-read.table(paste0("hepg2/data/star/",files_one_g[1]),header = FALSE,as.is = TRUE)[,1:7]
  
  
  for(i_file in files_one_g[-1]){
    
    files_all_one_gene<-read.table(paste0("hepg2/data/star/",i_file),header = FALSE,as.is = TRUE)[,1:7]
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
  colnames(sj_all)<-c("chr","start","end","strand","type","anno","uni_reads");
  
  sj_all<-sj_all[,c("chr","start","end","strand","type","uni_reads")];
    
  sj_all_mer<-sj_all %>% dplyr::group_by(chr,start,end,strand,type) %>% dplyr::summarise( anno=sum(uni_reads) );
  
  
  #x5_pos is the postion of 5'ss
  is_neg<-sj_all_mer[,4]==2;
  
  a<-sj_all_mer[is_neg,2];
  sj_all_mer[is_neg,2]<-sj_all_mer[is_neg,3];
  sj_all_mer[is_neg,3]<-a;
  
  
  colnames(sj_all_mer)<-c("chr","X5_pos","X3_pos","strand","type","anno");
  
  ###motif of 5'ss and 3'ss anchor
  sj_all_mer[(sj_all_mer$type==1),"X5_n"]<-"GT"
  sj_all_mer[(sj_all_mer$type==1),"X3_n"]<-"AG"
  
  sj_all_mer[(sj_all_mer$type==2),"X5_n"]<-"GT"
  sj_all_mer[(sj_all_mer$type==2),"X3_n"]<-"AG"
  
  sj_all_mer[(sj_all_mer$type==3),"X5_n"]<-"GC"
  sj_all_mer[(sj_all_mer$type==3),"X3_n"]<-"AG"
  
  sj_all_mer[(sj_all_mer$type==4),"X5_n"]<-"GC"
  sj_all_mer[(sj_all_mer$type==4),"X3_n"]<-"AG"
  
  sj_all_mer[(sj_all_mer$type==5),"X5_n"]<-"AT"
  sj_all_mer[(sj_all_mer$type==5),"X3_n"]<-"AC"
  
  sj_all_mer[(sj_all_mer$type==6),"X5_n"]<-"AT"
  sj_all_mer[(sj_all_mer$type==6),"X3_n"]<-"AC"
  
 
  sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno")];
  
  sj_all_mer_filter<-sj_all_mer;
  #sj_all_mer_filter<-sj_all_mer[sj_all_mer$X5_n=="GT",];
  #sj_all_mer_filter<-sj_all_mer[sj_all_mer$X3_n=="AG",];
  
  sj_all_mer_filter_sense<-sj_all_mer_filter[sj_all_mer_filter$X5_pos!=sj_all_mer_filter$X3_pos,];
  
  #g<-str_c("rRNA_",g);
  
  write.table(sj_all_mer_filter_sense,file = paste0("hepg2/data/star_mer/",g,"_tab"),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = TRUE );
  
}


################################################merge all ctl############################################

#   "ENCSR109IQO_CTL_RNA"

files_all<-list.files("/Users/mengli/Documents/projects/abs/hepg2/data/star_mer/");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);
gene_ids_ctl<-gene_ids[str_detect(gene_ids,"CTL_")];

#gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_norm_4");
#""
flag_id<-"aaa"
gene_ids_ctl<-setdiff(gene_ids_ctl,flag_id);

gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_all");

sj_all<-read.table(paste0("hepg2/data/star_mer/",gene_ids_ctl[1],"_tab"),header = TRUE,as.is = TRUE,sep = "\t");

sj_all<-cbind(sj_all,rep(str_sub(gene_ids_ctl[1],1,11),nrow(sj_all))  );

colnames(sj_all)[ncol(sj_all)]<-"rbp_count";

for( g in gene_ids_ctl[-1]){
  print(g)
  sj_t<-read.table(paste0("hepg2/data/star_mer/",g,"_tab"),header = TRUE,as.is = TRUE,sep = "\t");
  sj_t<-cbind(sj_t,rep(str_sub(g,1,11),nrow(sj_t) ) ) ;
  
  colnames(sj_t)[ncol(sj_t)]<-"rbp_count"
  
  sj_all<-rbind(sj_all,sj_t);
  
}


#sj_all_mer<-ddply(sj_all,.(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n),function(x){
#  c(paste0(x[,"anno"],collapse = " "), length(unique(x[,"rbp_count"]) )  );
#});

#colnames(sj_all_mer)[(ncol(sj_all_mer)-1):ncol(sj_all_mer)]<-c("anno","rbp_count");

#sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
#  summarise(anno=paste0(anno,collapse = " "),rbp_count=n_distinct(rbp_count ));
sj_all$rbp_count<-as.factor(sj_all$rbp_count);

sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
  dplyr::summarise(anno=sum(anno),rbp_count=n());

#sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
#  dplyr::summarise(anno=sum(anno));

sj_all_mer<-as.data.frame(sj_all_mer);
sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno","rbp_count")];

write.table(sj_all_mer,file = paste0("hepg2/data/star_mer/CTL_all_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );


#sj_all_mer_k562<-read.table(file = paste0("data/star_mer/CTL_all_tab"),header = TRUE,sep = "\t",as.is = TRUE)
    
#sj_all_mer<-rbind(sj_all_mer,sj_all_mer_k562);      

#write.table(sj_all_mer,file = paste0("hepg2/data/star_mer/CTL_all_tab2"),quote = FALSE,sep = "\t",
#            row.names = FALSE,col.names = TRUE );



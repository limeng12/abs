setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);

#   find /picb/rnasys2/limeng/fly/star_fly/ -type f -name '*SJ.out.tab' -exec cp '{}' /picb/rnasys2/limeng/fly/star_fly_re/ ';'
#   find /picb/rnasys2/limeng/fly/star_fly/ -type f -name '*.final.out' -exec cp '{}' /picb/rnasys2/limeng/fly/star_fly_log/ ';'

#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/fly/star_fly_re/\* /Users/mengli/Documents/projects/abs/fly/data/star_fly/
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/fly/star_fly_log/\* /Users/mengli/Documents/projects/abs/fly/data/star_log_fly/

#gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt"),readLines("fly/samples/gene_id_fly_all.txt") ) ));

gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt") ) ));


flag_id<-"nothing";

files_all<-list.files("fly/data/star_fly/");

#files<-sapply(str_split(files_all,"\\."),"[",1);

#gene_ids<-"normal"
##min overhang 8

for( g in gene_ids){
  print(g);
  
  
  files_one_g<-files_all[str_detect(files_all,g)];
  
  #sj_all<-read.table(paste0("fly/data/star_fly/",files_one_g[1]),header = FALSE,as.is = TRUE)
  
  #c("chr","start","end","strand","type","anno","uni_reads","mul_reads","overhang");
  sj_all<-data.frame(chr=c(),start=c(),end=c(),strand=c(),type=c(),anno=c(),uni_reads=c(),mul_reads=c(),overhang=c())
    
  for(i_file in files_one_g){
    
    if(!file.exists(paste0("fly/data/star_fly/",i_file)) | file.size(paste0("fly/data/star_fly/",i_file) )==0 ) {
      print(paste0("don't have ",i_file) );
      next;
    }
    
    
    files_all_one_gene<-read.table(paste0("fly/data/star_fly/",i_file),header = FALSE,as.is = TRUE)
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
  colnames(sj_all)<-c("chr","start","end","strand","type","anno","uni_reads","mul_reads","overhang");
  
  sj_all_mer<-sj_all %>% dplyr::group_by(chr,start,end,strand,type) %>% dplyr::summarise(anno=sum(uni_reads) );
  
  
  
  #x5_pos is the postion of 5'ss
  is_neg<-sj_all_mer[,4]==2;
  
  a<-sj_all_mer[is_neg,2];
  sj_all_mer[is_neg,2]<-sj_all_mer[is_neg,3];
  sj_all_mer[is_neg,3]<-a;
  
  #for( i in 1:nrow(sj_all_mer) ){
  #  if(sj_all_mer[i,4]==2){
  #    a<-sj_all_mer[i,2]
  #    sj_all_mer[i,2]<-sj_all_mer[i,3]
  #    sj_all_mer[i,3]<-a
  #  }
    
  #}
  
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
  
  write.table(sj_all_mer_filter_sense,file = paste0("fly/data/star_mer/",g,"_tab"),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = TRUE );
  
}


################################################merge all ctl############################################

#files_all<-list.files("/Users/mengli/Documents/projects/abs/data/star_mer/");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
#gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);


gene_ids_ctl<-gene_ids[str_detect(gene_ids,"Control") | str_detect(gene_ids,"CTL")];

#gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_norm_4");
#""
gene_ids_ctl<-setdiff(gene_ids_ctl,flag_id);

gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_all");

sj_all<-read.table(paste0("fly/data/star_mer/",gene_ids_ctl[1],"_tab"),header = TRUE,as.is = TRUE,sep = "\t");

sj_all<-cbind(sj_all,rep(gene_ids_ctl[1],nrow(sj_all))  );

colnames(sj_all)[ncol(sj_all)]<-"rbp_count";

for( g in gene_ids_ctl[-1]){
  print(g)
  sj_t<-read.table(paste0("fly/data/star_mer/",g,"_tab"),header = TRUE,as.is = TRUE,sep = "\t");
  sj_t<-cbind(sj_t,rep(g,nrow(sj_t) ) ) ;
  
  colnames(sj_t)[ncol(sj_t)]<-"rbp_count"
  
  sj_all<-rbind(sj_all,sj_t);
  
}


sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
  dplyr::summarise(anno=paste0(anno,collapse = " "),rbp_count=n() );


colnames(sj_all_mer)[(ncol(sj_all_mer)-1):ncol(sj_all_mer)]<-c("anno","rbp_count");

sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno","rbp_count")];

write.table(sj_all_mer,file = paste0("fly/data/star_mer/CTL_all_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );



# sj_all_mer<-ddply(sj_all,.(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n),function(x){
#   
#   c(paste0(x[,"anno"],collapse = " "), length(unique(x[,"rbp_count"]) )  );
#   
# });

#sj_all_mer<-sj_all %>% group_by(chr,start,end,strand,type) %>% summarise(anno=sum(uni_reads));
#sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
#  dplyr::summarise(anno=paste0(anno,collapse = " "),rbp_count=length(unique(rbp_count) ));



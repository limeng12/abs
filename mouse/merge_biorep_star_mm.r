setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);


#   find star_mm/ -type f -name '*SJ.out.tab' -exec cp '{}' star_re/ ';'
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/mouse/star_re/\* /Users/mengli/Documents/projects/abs/mouse/star/

#   find star_mm/ -type f -name '*.final.out' -exec cp '{}' star_log/ ';'
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/mouse/star_log/\* /Users/mengli/Documents/projects/abs/mouse/star_log/

gene_ids<-(unique(readLines("mouse/gene_id_mm.txt") ) );
#gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

#gene_ids<-c(gene_ids_bren,gene_ids_sicr);

#system("cp data/star/* data/star_sicr/");
#gene_ids<-unique(readLines("samples/gene_id_sicr.txt") );

flag_id<-"nothing"


#gene_ids<-"normal"
##min overhang 8

files_all<-list.files("mouse/star/");

#gene_ids<-"ENCSR443QFD_CTL_RNA_RNAi_inte"
for( g in gene_ids){
  print(g);
  
  files_one_g<-files_all[str_detect(files_all,fixed(g) )];
  
  if(length(files_one_g)==0){
    print(str_c("fail: ",g))
    next;
  }
  
  sj_all<-read.table(paste0("mouse/star/",files_one_g[1]),header = FALSE,as.is = TRUE)[,1:7]
  
  
  for(i_file in files_one_g[-1]){
    
    files_all_one_gene<-read.table(paste0("mouse/star/",i_file),header = FALSE,as.is = TRUE)[,1:7]
    
    sj_all<-rbind(sj_all,files_all_one_gene);
  }
  
  colnames(sj_all)<-c("chr","start","end","strand","type","anno","uni_reads");
  
  sj_all<-sj_all[,c("chr","start","end","strand","type","uni_reads")];
  
  
  #colnames(sj_all)<-c("chr","start","end","strand","type","anno","uni_reads","mul_reads","overhang");
  #sj_all[,"sample_count"]<-length(files_one_g);
  #sj_all_mer<-ddply(sj_all,.(chr,start,end,strand,type),.fun=function(x){
    #type_s<-paste(x[,"type"],collapse = ",");
  #  uni_reads_c<-sum(x[,"uni_reads"]);
  #  uni_reads_c
    #paste0("type=",type_s," unique_read_c=",uni_reads_c);
  #},.parallel=TRUE );
  #sj_all_mer<-ddply(sj_all,.(chr,start,end,strand,type),summarize, uni_reads_c=sum(uni_reads) );
    
  sj_all_mer<-sj_all %>% group_by(chr,start,end,strand,type) %>% summarise(anno=sum(uni_reads),sample_count=length(uni_reads) );
  
  #minimium number reads for control is 1, for rbp is 1
  if(str_detect(g,"CTL_") && (g!=flag_id)  ) {
    sj_all_mer<-sj_all_mer[sj_all_mer$anno>0,];
  }else{
    sj_all_mer<-sj_all_mer[sj_all_mer$anno>0,];
  }
  
  
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
  
  colnames(sj_all_mer)<-c("chr","X5_pos","X3_pos","strand","type","anno","sample_count");
  
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
  
 
  sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno","sample_count")];
  
  sj_all_mer_filter<-sj_all_mer;
  #sj_all_mer_filter<-sj_all_mer[sj_all_mer$X5_n=="GT",];
  #sj_all_mer_filter<-sj_all_mer[sj_all_mer$X3_n=="AG",];
  
  sj_all_mer_filter_sense<-sj_all_mer_filter[sj_all_mer_filter$X5_pos!=sj_all_mer_filter$X3_pos,];
  
  #sj_all_mer_filter_sense[,"sample_count"]<-length(files_one_g);
  #g<-str_c("rRNA_",g);
  
  write.table(sj_all_mer_filter_sense,file = paste0("mouse/star_mer/",g,"_tab"),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = TRUE );
  
}


################################################merge all ctl############################################

#files_all<-list.files("mouse/star_mer/");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
#gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);
#gene_ids<-(unique(readLines("mouse/gene_id_mm.txt") ) );


gene_ids<-(unique(readLines("mouse/gene_id_mm.txt") ) );

#files_all<-list.files("mouse/star_mer/");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
#gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);


gene_ids_ctl<-gene_ids[str_detect(gene_ids,"siNT") | str_detect(gene_ids,"untreated") | str_detect(gene_ids,"mock")];

gene_ids_ctl<-gene_ids_ctl[str_detect(gene_ids_ctl,"CGR8")];
#gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_norm_4");
#""
gene_ids_ctl<-setdiff(gene_ids_ctl,flag_id);

gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_all");

sj_all<-read.table(paste0("mouse/star_mer/",gene_ids_ctl[1],"_tab"),header = TRUE,as.is = TRUE,sep = "\t");

#sj_all<-cbind(sj_all,rep(gene_ids_ctl[1],nrow(sj_all))  );

#colnames(sj_all)[ncol(sj_all)]<-"rbp_count";

for( g in gene_ids_ctl[-1]){
  print(g)
  sj_t<-read.table(paste0("mouse/star_mer/",g,"_tab"),header = TRUE,as.is = TRUE,sep = "\t");
  #sj_t<-cbind(sj_t,rep(g,nrow(sj_t) ) ) ;
  
  #colnames(sj_t)[ncol(sj_t)]<-"rbp_count"
  
  sj_all<-rbind(sj_all,sj_t);
  
}

#sj_all_mer<-sj_all %>% group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
#  summarise(anno=paste0(anno,collapse = " "),rbp_count=length(unique(sample_count) ));

sj_all_mer<-sj_all %>% dplyr::group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
  dplyr::summarise(anno=paste0(anno,collapse = " "),sample_count=sum(sample_count) );


sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno","sample_count")];

write.table(sj_all_mer,file = paste0("mouse/star_mer/CTL_all_CGR8_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );




gene_ids<-(unique(readLines("mouse/gene_id_mm.txt") ) );

#files_all<-list.files("mouse/star_mer/");
#cut -f1-3,5 star_target_only_jc/AQR_no_ctl.bed > star_target_only_jc_bedgraph/AQR_no_ctl.bedgraph
#gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);


gene_ids_ctl<-gene_ids[str_detect(gene_ids,"siNT") | str_detect(gene_ids,"untreated") | str_detect(gene_ids,"mock")];

gene_ids_ctl<-gene_ids_ctl[str_detect(gene_ids_ctl,"N2A")];
#gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_norm_4");
#""
gene_ids_ctl<-setdiff(gene_ids_ctl,flag_id);

gene_ids_ctl<-setdiff(gene_ids_ctl,"CTL_all");

sj_all<-read.table(paste0("mouse/star_mer/",gene_ids_ctl[1],"_tab"),header = TRUE,as.is = TRUE,sep = "\t");

#sj_all<-cbind(sj_all,rep(gene_ids_ctl[1],nrow(sj_all))  );

#colnames(sj_all)[ncol(sj_all)]<-"rbp_count";

for( g in gene_ids_ctl[-1]){
  print(g)
  sj_t<-read.table(paste0("mouse/star_mer/",g,"_tab"),header = TRUE,as.is = TRUE,sep = "\t");
  #sj_t<-cbind(sj_t,rep(g,nrow(sj_t) ) ) ;
  
  #colnames(sj_t)[ncol(sj_t)]<-"rbp_count"
  
  sj_all<-rbind(sj_all,sj_t);
  
}

#sj_all_mer<-sj_all %>% group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
#  summarise(anno=paste0(anno,collapse = " "),rbp_count=length(unique(rbp_count) ));

sj_all_mer<-sj_all %>% group_by(chr,X5_pos,X3_pos,strand,type,X5_n,X3_n) %>% 
  summarise(anno=paste0(anno,collapse = " "),sample_count=sum(sample_count) );



sj_all_mer<-sj_all_mer[,c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno","sample_count")];

write.table(sj_all_mer,file = paste0("mouse/star_mer/CTL_all_N2A_tab"),quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = TRUE );




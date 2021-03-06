setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);
library(BSgenome);
library(BSgenome.Dmelanogaster.UCSC.dm6);
library(ggplot2);
library(RWebLogo);
options("scipen"=100, "digits"=4);

dm_genome <- getBSgenome("BSgenome.Dmelanogaster.UCSC.dm6");

#gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt"),readLines("fly/samples/gene_id_fly_all.txt") ) ));

gene_ids<-sort(unique(c(readLines("fly/samples/gene_id_fly.txt") ) ));


ctl_sj<-read.table("fly/data/star_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);

#colnames(ctl_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno");

###high confidence control should have enought unique read mapping and also in annotation.
##for control calculation only
ctl_sj_high_confi<-ctl_sj[ctl_sj$anno>2,];


## delete all introns' 5'ss and 3'ss in ENCODE
anno_sj<-read.table("fly/anno/dm6_ensembl_intron_from_ucsc.bed",sep = "\t",header = FALSE, as.is = TRUE);
#anno_sj<-anno_sj[!duplicated(anno_sj[,c("chrom","start","end")]  ),];
anno_sj[,2]<-anno_sj[,2]+1;
anno_sj[,1]<-str_sub(anno_sj[,1],4);


anno_sj<-unname(anno_sj);

anno_sj_positive<-anno_sj[anno_sj[,6]=="+",];
anno_sj_negative<-anno_sj[anno_sj[,6]=="-",];

anno_5ss<-unique( c(str_c(anno_sj_positive[,1],anno_sj_positive[,2]),
            str_c(anno_sj_negative[,1],anno_sj_negative[,3])  )  );


anno_3ss<-unique( c(str_c(anno_sj_positive[,1],anno_sj_positive[,3]),
            str_c(anno_sj_negative[,1],anno_sj_negative[,2])  ) );


#total_5ss<-unique(intersect(str_c(ctl_sj_high_confi[,"chr"],ctl_sj_high_confi[,"X5_pos"] ), anno_5ss) );
#total_3ss<-unique(intersect(str_c(ctl_sj_high_confi[,"chr"],ctl_sj_high_confi[,"X3_pos"] ), anno_3ss) );

total_5ss<-anno_5ss;

total_3ss<-anno_3ss;


gene_ids_noctl<-gene_ids[!str_detect(gene_ids,"CTL_all")];

junc_len<-matrix(nrow=0,ncol=8);

g_53_sj<-matrix(nrow=0,ncol=9);



for(g in gene_ids){
  print(g);
  if(str_detect(g,"_CTL")){
    next;
  }
  
  
  t_sj<-read.table(paste0("fly/data/star_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  
  colnames(t_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n","anno");
  
  t_sj<-t_sj[t_sj$anno>1,];
  
  #t_sj<-t_sj[str_detect(t_sj$chr,"chr"),];
  
  
  t_sj$strand2<-sapply(t_sj$strand,function(x){
    if(x==1){
      return("+");
    }
      return("-");
  });
  
  t_sj_5only_nondup_index<-!duplicated(t_sj[,c("chr","X5_pos","X5_n")] );
  
  
  #t_5_only_sj_count<-length(setdiff(t_sj[,2],ctl_sj[,2]) );
  t_5_only_sj_set_index<-(str_c(t_sj[,"chr"],t_sj[,"X5_pos"],t_sj[,"X5_n"]) %in%
                           str_c(ctl_sj[,"chr"],ctl_sj[,"X5_pos"],ctl_sj[,"X5_n"])  ) |
                             (str_c(t_sj[,"chr"],t_sj[,"X5_pos"] ) %in%
                                anno_5ss   );
  #t_5_only_sj_count<-sum(!t_5_only_sj_set_index);
  
  ctl_5_only_sj_count<-length(unique(setdiff(total_5ss, str_c(t_sj[,"chr"],t_sj[,"X5_pos"])  )  ) );
  
  #t_5_cano_per<-sum(str_detect(t_5_only_sj_set,"GT")) /length(t_5_only_sj_set);
  t_5_only_sj_set<-unique(str_c(t_sj[,"chr"],t_sj[,"X5_pos"],t_sj[,"X5_n"])[!t_5_only_sj_set_index]);
  t_5_only_sj_count_u12<-length(t_5_only_sj_set[str_detect(t_5_only_sj_set,"AT")] );
  
  
  t_5_only_sj_count_gc<-length(t_5_only_sj_set[str_detect(t_5_only_sj_set,"GC")] );
  t_5_only_sj_count<-length(t_5_only_sj_set[str_detect(t_5_only_sj_set,"GT")] );
  
  t_sj_in_only_5ss_fr<-t_sj[!t_5_only_sj_set_index,c("chr","X5_pos","X5_pos","strand")];
  
  if(nrow(t_sj_in_only_5ss_fr)>0){
    t_sj_in_only_5ss_fr[,2]<-t_sj_in_only_5ss_fr[,2]-1;
    
    write.table(cbind(t_sj_in_only_5ss_fr,1),file=paste0("fly/data/star_target_only_jc_5ss/",g,"_no_ctl.bed"), sep="\t",
                quote = FALSE,col.names = FALSE,row.names = FALSE);
  }
  
  positive_index<-(t_sj$strand2=="+" ) & 
    (str_detect(t_sj$chr,fixed("2L"))|
       str_detect(t_sj$chr,fixed("2R")) | str_detect(t_sj$chr,fixed("3L"))  | str_detect(t_sj$chr,fixed("3R")) |
        str_detect(t_sj$chr,fixed("Y")) | str_detect(t_sj$chr,fixed("X")) ) & (!str_detect(t_sj$chr,"Y_"));
  negative_index<-(t_sj$strand2=="-" ) & 
    (str_detect(t_sj$chr,fixed("2L"))|
       str_detect(t_sj$chr,fixed("2R")) | str_detect(t_sj$chr,fixed("3L"))  | str_detect(t_sj$chr,fixed("3R")) |
       str_detect(t_sj$chr,fixed("Y")) | str_detect(t_sj$chr,fixed("X")) ) & (!str_detect(t_sj$chr,"Y_"));
  
  
  if(sum(!t_5_only_sj_set_index)>0){
     
     t_5_seqs<-getSeq(dm_genome,str_c("chr",t_sj[(!t_5_only_sj_set_index)&positive_index,"chr"]),
                      t_sj[(!t_5_only_sj_set_index)&positive_index,"X5_pos"]-3,
                      t_sj[(!t_5_only_sj_set_index)&positive_index,"X5_pos"]+5,
                      strand=t_sj[(!t_5_only_sj_set_index)&positive_index,"strand2"],as.character=TRUE);
     #cat(as.character(t_5_seqs),file=paste0("fly/data/star_abs5_motif/",g,".seq"), sep="\n");
     
     t_5_seqs_neg<-getSeq(dm_genome,str_c("chr",t_sj[(!t_5_only_sj_set_index)&(negative_index),"chr"]),
                      t_sj[(!t_5_only_sj_set_index)&(negative_index),"X5_pos"]-5,
                      t_sj[(!t_5_only_sj_set_index)&(negative_index),"X5_pos"]+3,
                      strand=t_sj[(!t_5_only_sj_set_index)&(negative_index),"strand2"],as.character=TRUE);
     #cat(as.character(t_5_seqs_neg),file=paste0("fly/data/star_abs5_motif/",g,".seq"), sep="\n",append=TRUE);
     
     t_5_seqs<-c(t_5_seqs,t_5_seqs_neg)
     
     weblogo(as.character(t_5_seqs),file.out=paste0("fly/data/star_abs5_motif/",g,".pdf") ,format="pdf" ,open = FALSE )
     
   }
  
  
  t_sj_3only_nondup_index<-!duplicated(t_sj[,c("chr","X3_pos","X3_n")] );
  
  #t_3_only_sj_count<-length(setdiff(t_sj[,3],ctl_sj[,3]) );
  t_3_only_sj_set_index<-(str_c(t_sj[,"chr"],t_sj[,"X3_pos"],t_sj[,"X3_n"]) %in%
                           str_c(ctl_sj[,"chr"],ctl_sj[,"X3_pos"],ctl_sj[,"X3_n"]) ) |
                              (str_c(t_sj[,"chr"],t_sj[,"X3_pos"] ) %in%
                                 anno_3ss  )
  #t_3_only_sj_count<-sum(!t_3_only_sj_set_index);
  
  t_3_only_sj_set<-unique(str_c(t_sj[,"chr"],t_sj[,"X3_pos"],t_sj[,"X3_n"])[!t_3_only_sj_set_index])
                           
  t_3_only_sj_count_u12<-length(t_3_only_sj_set[str_detect(t_3_only_sj_set,"AC")] );
  
  #t_3_only_sj_count_u12_readcount<-sum(t_3_only_sj_set[str_detect(t_3_only_sj_set,"AC")])
    
  
  t_3_only_sj_count<-length(t_3_only_sj_set[str_detect(t_3_only_sj_set,"AG")] );
  
  positive_index<-(t_sj$strand2=="+")  & 
    (str_detect(t_sj$chr,fixed("2L")) |
       str_detect(t_sj$chr,fixed("2R")) | str_detect(t_sj$chr,fixed("3L"))  | str_detect(t_sj$chr,fixed("3R")) |
       str_detect(t_sj$chr,fixed("Y")) | str_detect(t_sj$chr,fixed("X")) ) & (!str_detect(t_sj$chr,"Y_"));
  
  negative_index<-(t_sj$strand2=="-")  & 
    (str_detect(t_sj$chr,fixed("2L")) |
       str_detect(t_sj$chr,fixed("2R")) | str_detect(t_sj$chr,fixed("3L"))  | str_detect(t_sj$chr,fixed("3R")) |
       str_detect(t_sj$chr,fixed("Y")) | str_detect(t_sj$chr,fixed("X")) ) & (!str_detect(t_sj$chr,"Y_"));
  
  
  if(sum(!t_3_only_sj_set_index)>0){
    t_3_seqs<-getSeq(dm_genome,str_c("chr",t_sj[(!t_3_only_sj_set_index)&positive_index,"chr"]),
                      t_sj[(!t_3_only_sj_set_index)&positive_index,"X3_pos"]-19,
                      t_sj[(!t_3_only_sj_set_index)&positive_index,"X3_pos"]+3,
                      strand=t_sj[(!t_3_only_sj_set_index)&positive_index,"strand2"],as.character=TRUE);
    #cat(as.character(t_3_seqs),file=paste0("fly/data/star_abs3_motif/",g,".seq"), sep="\n");
    
    
    t_3_seqs_neg<-getSeq(dm_genome,str_c("chr",t_sj[(!t_3_only_sj_set_index)&(negative_index),"chr"]),
                     t_sj[(!t_3_only_sj_set_index)&(negative_index),"X3_pos"]-3,
                     t_sj[(!t_3_only_sj_set_index)&(negative_index),"X3_pos"]+19,
                     strand=t_sj[(!t_3_only_sj_set_index)&(negative_index),"strand2"],as.character=TRUE);
    #cat(as.character(t_3_seqs_neg),file=paste0("fly/data/star_abs3_motif/",g,".seq"), sep="\n",append = TRUE);
    
    t_3_seqs<-c(t_3_seqs,t_3_seqs_neg);
    

    weblogo(as.character(t_3_seqs),file.out=paste0("fly/data/star_abs3_motif/",g,".pdf") ,format="pdf" ,open = FALSE )
    
   }
   
  
  t_sj_in_only_3ss_fr<-t_sj[!t_3_only_sj_set_index,c("chr","X3_pos","X3_pos","strand"),];
  
  if(nrow(t_sj_in_only_3ss_fr)>0){
  
    t_sj_in_only_3ss_fr[,2]<-t_sj_in_only_3ss_fr[,2]-1;
  
    write.table(cbind(t_sj_in_only_3ss_fr,1),file=paste0("fly/data/star_target_only_jc_3ss/",g,"_no_ctl.bed"), sep="\t",
              quote = FALSE,col.names = FALSE,row.names = FALSE);
    
  }
  
  #t_5_only_sj_set_index_all<-(str_c(t_sj[,"chr"],t_sj[,"X5_pos"],t_sj[,"X5_n"]) %in%
  #                          str_c(ctl_sj[,"chr"],ctl_sj[,"X5_pos"],ctl_sj[,"X5_n"]) );
  #t_3_only_sj_set_index_all<-(str_c(t_sj[,"chr"],t_sj[,"X3_pos"],t_sj[,"X3_n"]) %in%
  #                          str_c(ctl_sj[,"chr"],ctl_sj[,"X3_pos"],ctl_sj[,"X3_n"]) );
  
  both_unique<-sum((!t_5_only_sj_set_index) & (!t_3_only_sj_set_index) );
  
  ctl_3_only_sj_count<-length(unique(setdiff(total_3ss,str_c(t_sj[,"chr"],t_sj[,"X3_pos"]) )  ));
  
  #t_3_cano_per<-sum(str_detect(t_3_only_sj_set,"AG"))/length(t_3_only_sj_set);
  
  
  
  t_ctl_sj<-left_join(t_sj,ctl_sj,
                      by=c("chr"="chr","X5_pos"="X5_pos","X3_pos"="X3_pos",
                           "strand"="strand","type"="type","X5_pos"="X5_pos","X3_pos"="X3_pos") );
  
  t_only_sj_all<-t_ctl_sj[is.na(t_ctl_sj$anno.y),];
  
  #junc_len[[g]]<-abs(t_only_sj_all$X3_pos-t_only_sj_all$X5_pos);
  
  junc_len<-rbind(junc_len,c(g,nrow(t_only_sj_all),
                             as.vector(summary(abs(t_only_sj_all$X3_pos-t_only_sj_all$X5_pos))) ) );
  
  
  is_neg<-t_only_sj_all[,"strand"]==2;
  
  a<-t_only_sj_all[is_neg,"X5_pos"];
  t_only_sj_all[is_neg,"X5_pos"]<-t_only_sj_all[is_neg,"X3_pos"];
  t_only_sj_all[is_neg,"X3_pos"]<-a;
  
  t_only_sj_all$X5_pos<-t_only_sj_all$X5_pos-1;
  
  if(nrow(t_only_sj_all)>0){
    t_only_sj_all$score<-1;
  }
  
  if("chrM" %in% t_only_sj_all$chr){
    write.table(t_only_sj_all[t_only_sj_all$chr=="chrM",],
                file=paste0("fly/data/star_wrong_jc/",g,"_no_ctl_somethingwrong.bed"), sep="\t",
                quote = FALSE,col.names = FALSE,row.names = FALSE);
  }
  
  
  write.table(t_only_sj_all,file=paste0("fly/data/star_target_only_jc/",g,"_no_ctl.bed"), sep="\t",
              quote = FALSE,col.names = FALSE,row.names = FALSE);
  
  
  #t_5_only_anno_sj_count<-length(setdiff(t_sj[,2],anno_sj[,2]) );
  #t_3_only_anno_sj_count<-length(setdiff(t_sj[,3],anno_sj[,3]) );
  #ctl_5_only_sj_count<-length(setdiff(ctl_sj[,2],t_sj[,2]) );
  #ctl_3_only_sj_count<-length(setdiff(ctl_sj[,3],t_sj[,3]) );
  
  g_53_sj<-rbind(g_53_sj,c(g,t_5_only_sj_count,t_3_only_sj_count,both_unique,
                           t_5_only_sj_count_u12,t_3_only_sj_count_u12,t_5_only_sj_count_gc,
                           ctl_5_only_sj_count,ctl_3_only_sj_count ) );
  #t_5_cano_per,t_3_cano_per,) )
  
}


g_53_sj_data<-as.data.frame(g_53_sj,stringsAsFactors=FALSE);
colnames(g_53_sj_data)<-c("rbp","cryptic 5'ss","cryptic 3'ss",
                          "cryptic whole junction","U12 cryptic 5'ss","U12 cryptic 3'ss",
                          "CG type 5'ss","ctrl_only_5ss","ctrl_only_3ss")
#,"X5_cano_per","X3_cano_per","X5_ctl","X3_ctl");


g_53_sj_data<-g_53_sj_data[order(g_53_sj_data[,"cryptic 3'ss"],decreasing = TRUE),];

write.table(g_53_sj_data,file="fly/result/g_53_sj_all_full_fly.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



gene_read_fr_fr<-read.table("fly/result/rbp_read_count_map.tsv",sep = "\t",as.is = TRUE,header = TRUE);


g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr_fr,by=c("rbp"="g") );

g_53_sj_data_read[,"corrected cryptic 3'ss"]<-as.numeric(g_53_sj_data_read[,"cryptic 3'ss"])/
  as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);

g_53_sj_data_read[,"corrected cryptic 5'ss"]<-as.numeric(g_53_sj_data_read[,"cryptic 5'ss"])/
  as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);

write.table(g_53_sj_data_read,file="fly/result/g_53_sj_all_read_count_fly.tsv",sep="\t",
             quote = FALSE,col.names = TRUE,row.names = FALSE);



#g_53_sj_data_read$X3_X5_ratio<-g_53_sj_data_read$X3/g_53_sj_data_read$X5;

#colnames(g_53_sj_data_read)<-c("rbp","cryptic 5'ss","cryptic 3'ss",
#                                "cryptic whole junction","U12 cryptic 5'ss","U12 cryptic 3'ss","CG type 5'ss",
#                                "read_count","cryptic 3'ss correct read count","cryptic 3'ss/cryptic 5'ss");


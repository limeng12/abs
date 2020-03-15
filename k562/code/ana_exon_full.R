setwd("/Users/mengli/Documents/projects/abs");
#library(plyr);
library(stringr);
library(dplyr);
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19);
options(scipen=999)

#/picb/rnasys/program/install/java/jdk1.8.0_45/bin/java -jar midec.jar exon/CTL_ENCSR000AEO_1_exon \
#star/CTL_ENCSR000AEO/CTL_ENCSR000AEO_1Aligned.sortedByCoord.out.bam
#/picb/rnasys/program/install/java/jdk1.8.0_45/bin/java -jar midec.jar exon/CTL_ENCSR000AEO_2_exon \
#star/CTL_ENCSR000AEO/CTL_ENCSR000AEO_2Aligned.sortedByCoord.out.bam

region_len<-function(x){
  region_list<-str_split(x,":|-");
  start<-as.numeric(sapply(region_list,"[",2));
  end<-as.numeric(sapply(region_list,"[",3));
  len<-end-start+1;
  return(len);
}

hg19_genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19");

#  getSeq(hg19_genome,data_filter[,"Region"],data_filter[,"Position"]-49,data_filter[,"Position"]+50);

#gene_ids<-(readLines("samples/gene_id.txt") );
gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

gene_ids<-c(gene_ids_bren,gene_ids_sicr);
gene_ids<-gene_ids[!str_detect(gene_ids,"_inte")]
#gene_ids<-setdiff(gene_ids,"ENCSR000KYM_FLAG_RNA");


anno_sj<-read.table("anno/hg19_gencode_intron_from_ucsc.bed",sep = "\t",header = FALSE, as.is = TRUE);
anno_sj[,2]<-anno_sj[,2]+1
anno_sj<-unname(anno_sj);

anno_sj_positive<-anno_sj[anno_sj[,6]=="+",];
anno_sj_negative<-anno_sj[anno_sj[,6]=="-",];


anno_5ss<-unique( c(str_c(anno_sj_positive[,1],":",anno_sj_positive[,2]),
                    str_c(anno_sj_negative[,1],":",anno_sj_negative[,3])  )  );


anno_3ss<-unique( c(str_c(anno_sj_positive[,1],anno_sj_positive[,3]),
                    str_c(anno_sj_negative[,1],anno_sj_negative[,2])  ) );


ctl_sj_sj_only<-read.table("data/star_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);

ctl_sj_5ss<-str_c(ctl_sj_sj_only$chr,":",ctl_sj_sj_only$X5_pos);

ctl_sj_3ss<-str_c(ctl_sj_sj_only$chr,":",ctl_sj_sj_only$X3_pos);



ctl_sj<-read.table("data/exon_mer/CTL_all_tab",sep = "\t",header = TRUE, as.is = TRUE);

ctl_sj_only_polyA_no_shCTL<-read.table("data/exon_mer/ENCSR000AEO_CTL_poly_tab",sep = "\t",header = TRUE, as.is = TRUE);


intron_coor<-read.table(file="anno/hg19_gencode_intron_from_ucsc.bed",sep="\t",
                        as.is = TRUE,header = FALSE);

intron_coor[,2]<-intron_coor[,2]+1;

colnames(intron_coor)<-c("chr","start","end","name","score","strand");

intron_coor_range<-with(intron_coor,
                        GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                                strand=strand) );


#exon_coor<-read.table(file="anno/hg19_refgene.gtf",sep="\t",as.is = FALSE,header = TRUE);
exon_coor<-read.table(file="anno/gencode.v29lift37.annotation_exon.gtf",sep="\t",as.is = FALSE,header = TRUE);

colnames(exon_coor)[1:7]<-c("chr","type","class","start","end","score","strand")
exon_coor<-unique(exon_coor[exon_coor[,3]=="exon",]);

all_exon_frame_0<-sum( (exon_coor[,"end"]-exon_coor[,"start"]+1 ) %% 3 == 0 );
all_exon_frame_1<-sum( (exon_coor[,"end"]-exon_coor[,"start"]+1 ) %% 3 == 1 );
all_exon_frame_2<-sum( (exon_coor[,"end"]-exon_coor[,"start"]+1 ) %% 3 == 2 );

print(paste0("all exon_frame = 0: ",all_exon_frame_0));
print(paste0("all exon_frame = 1: ",all_exon_frame_1));
print(paste0("all exon_frame = 2: ",all_exon_frame_2));


exon_coor_region<-str_c(exon_coor$chr,":",exon_coor$start,"-",exon_coor$end);
exon_coor_range<-with(exon_coor,
                      GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                              strand=strand) );



#anno_sj<-read.table("anno/intron_coor_gencode.tsv",sep = "\t",header = TRUE, as.is = TRUE);


gene_ids_noctl<-gene_ids[!str_detect(gene_ids,"CTL_all")];

junc_len<-matrix(nrow=0,ncol=8);

g_53_sj<-matrix(nrow=0,ncol=8);


for(g in gene_ids){
  print(g);
  
  if(!file.exists(paste0("data/exon_mer/",g,"_tab"))){
    print(paste0("don't have ",g))
    next;
  }
  
  gene_sj_sj_only<-read.table(paste0("data/star_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  
  gene_sj_3ss<-str_c(gene_sj_sj_only$chr,":",gene_sj_sj_only$X3_pos);
  
  
  
  t_sj<-read.table(paste0("data/exon_mer/",g,"_tab"),sep = "\t",header = TRUE, as.is = TRUE);
  
  ##at least two reads for one exon
  
  #t_sj<-t_sj[t_sj$read_count>1,];
  t_sj<-t_sj[!is.na(t_sj$region),];
  t_sj<-t_sj[!is.na(t_sj$read_count),];
  

  t_only_sj_all<-anti_join(t_sj,ctl_sj,
                      by=c("region"="region") ) #%>% anti_join(exon_coor_region,by=c("region"="region"));
  
  t_only_sj_all<-t_only_sj_all[!(t_only_sj_all$region %in% exon_coor_region),];
  
  
  all_cryptic_exon_len<-region_len(t_only_sj_all$region);
  
  number_of_0_frame<-sum(all_cryptic_exon_len %% 3 ==0);
  number_of_1_frame<-sum(all_cryptic_exon_len %% 3 ==1);
  number_of_2_frame<-sum(all_cryptic_exon_len %% 3 ==2);
  
  #ctl_sj_t<-ctl_sj[ctl_sj$read_count>5,];
  #ctl_sj_t<-ctl_sj[ctl_sj$gene>4,];
  #############
  ctl_sj_t<-ctl_sj_only_polyA_no_shCTL;
  #############
  t_only_ctl_sj<-anti_join(ctl_sj_t,t_sj, by=c("region"="region") );
  
  t_only_ctl_sj_small<-t_only_ctl_sj[region_len(t_only_ctl_sj$region)<50,];
  
  
  t_only_anno_ctl_sj<-anti_join(ctl_sj_t[ctl_sj_t$region %in% exon_coor_region,],t_sj,
                           by=c("region"="region") );
  
  t_only_anno_ctl_sj_small<-t_only_anno_ctl_sj[region_len(t_only_anno_ctl_sj$region)<50,]
  
  
  region_list<-str_split(t_only_sj_all$region,":|-");
  chr<-(sapply(region_list,"[",1));
  
  start<-as.numeric(sapply(region_list,"[",2) );
  end<-as.numeric(sapply(region_list,"[",3) );
  len<-end-start+1;
  
  
  t_only_sj_all_bed<-data.frame(chr=chr,
                                start=start,
                                end=end,
                                read_count=t_only_sj_all$read_count);
  
  if(nrow(t_only_sj_all_bed)>0){
    t_only_sj_all_bed<-t_only_sj_all_bed[setdiff(1:nrow(t_only_sj_all_bed),which(is.na(t_only_sj_all_bed[,1])) ), ]
    
    
    write.table(cbind(t_only_sj_all_bed,1),
                file=paste0("data/exon_mer_target_only/",g,"_no_ctl.bed"), sep="\t",
                quote = FALSE,col.names = FALSE,row.names = FALSE);
    
    
    t_only_sj_all_bed_range<-with(t_only_sj_all_bed,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                                                            strand="*",read_count=read_count) );
    
    
    t_only_sj_all_bed_range_count_intron<-countOverlaps(t_only_sj_all_bed_range,intron_coor_range,type="within");
    
    
    t_only_sj_all_bed_output<-t_only_sj_all_bed[t_only_sj_all_bed_range_count_intron>0,];
    
    t_only_sj_all_bed_output$start<-t_only_sj_all_bed_output$start-1;
    
    write.table(cbind(t_only_sj_all_bed_output,1),
                file=paste0("data/exon_mer_target_only_intron/",g,"_no_ctl.bed"), sep="\t",
                quote = FALSE,col.names = FALSE,row.names = FALSE);
    
    
    
    
    t_only_sj_all_bed_pos<-data.frame(chr=chr,
                                  start=end,
                                  end=end+1,
                                  read_count=t_only_sj_all$read_count,
                                  strand=rep("+",length(chr))  );
    
    t_only_sj_all_bed_neg<-data.frame(chr=chr,
                                      start=start-2,
                                      end=start-1,
                                      read_count=t_only_sj_all$read_count,
                                      strand=rep("-",length(chr)));
    
    t_only_sj_all_bed_all<-cbind(rbind(t_only_sj_all_bed_pos,t_only_sj_all_bed_neg),1);
    
    t_only_sj_all_bed_all[,"may_5ss"]<-str_c(t_only_sj_all_bed_all$chr,":",t_only_sj_all_bed_all$end);
    
    
    ###cryptic exon induced cryptic 5ss
    ##remove 5ss in control and annotation, also remove 3ss
    ##cryptic exons' 3ss must be in total cryptic 3ss
    ##5ss and 3ss must not be overlapped
    t_only_sj_all_bed_all_cryptic_5ss<-t_only_sj_all_bed_all[!(t_only_sj_all_bed_all$may_5ss %in% c(ctl_sj_5ss,ctl_sj_3ss,
                                                                                                    anno_5ss,anno_3ss,
                                                                                                    gene_sj_3ss) ),];
    
    write.table(t_only_sj_all_bed_all_cryptic_5ss,
                file=paste0("data/exon_mer_target_only_5ss/",g,"_no_ctl.bed"), sep="\t",
                quote = FALSE,col.names = FALSE,row.names = FALSE);
    
    
    g_53_sj<-rbind(g_53_sj,c(g,number_of_0_frame,number_of_1_frame,number_of_2_frame,
                             #sum(t_only_sj_all_bed_small_range_count_intron>0),
                             #nrow(t_only_sj_all_bed_small),
                             nrow(t_only_sj_all),
                             #####nrow(ctl_sj)-
                               nrow(t_only_ctl_sj),
                             #nrow(t_only_ctl_sj_small),
                             
                             #####length(exon_coor_region)-
                             nrow(t_only_anno_ctl_sj),
                             #nrow(t_only_anno_ctl_sj_small),
                             sum(t_only_sj_all$read_count,na.rm = TRUE) ) );

  }else{
    g_53_sj<-rbind(g_53_sj,c(g,number_of_0_frame,number_of_1_frame,number_of_2_frame,0,0,0,0) );
  }
  
}


g_53_sj_data<-data.frame(rbp=g_53_sj[,1],
                         number_of_0_frame=g_53_sj[,2],
                         number_of_1_frame=g_53_sj[,3],
                         number_of_2_frame=g_53_sj[,4],
                         #no_in_ctl_exon_count_small_exon=g_53_sj[,5],
                         #no_in_ctl_exon_count_small_intron=g_53_sj[,6],
                         #no_in_ctl_exon_count_small=g_53_sj[,7],
                         no_in_ctl_exon_count=g_53_sj[,5],
                         missed_exon_insh_count=g_53_sj[,6],
                         #missed_exon_small_insh_count=g_53_sj[,10],
                         
                         missed_anno_exon_insh_count=g_53_sj[,7],
                         #missed_anno_small_exon_insh_count=g_53_sj[,12],
                         
                         total_read_count_support=g_53_sj[,8]);


g_53_sj_data<-unique(g_53_sj_data);

g_53_sj_data<-g_53_sj_data[order(g_53_sj_data$no_in_ctl_exon_count,decreasing = TRUE),];

write.table(g_53_sj_data,file="result/g_53_exon_all_full.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);


gene_read_fr_fr<-read.table("result/rbp_read_count_map.tsv",sep = "\t",as.is = TRUE,header = TRUE);

g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr_fr,by=c("rbp"="g") );


g_53_sj_data_read[,"corrected missed_exon_insh_count"]<-
  as.numeric(as.character(g_53_sj_data_read[,"missed_exon_insh_count"]))/
  as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);

write.table(g_53_sj_data_read,file="result/g_53_exon_all_full_only_read_count.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);



# g_53_sj_data_read[,"corrected no_in_ctl_exon_count"]<-
#   as.numeric(as.character(g_53_sj_data_read[,"no_in_ctl_exon_count"]))/
#   as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);

# g_53_sj_data_read[,"corrected no_in_ctl_exon_count_small"]<-
#   as.numeric(as.character(g_53_sj_data_read[,"no_in_ctl_exon_count_small"]))/
#   as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);
# 
# g_53_sj_data_read[,"corrected no_in_ctl_exon_count_small_intron"]<-
#   as.numeric(as.character(g_53_sj_data_read[,"no_in_ctl_exon_count_small_intron"]))/
#   as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);

# 
# g_53_sj_data_read[,"corrected missed_exon_small_insh_count"]<-
#   as.numeric(as.character(g_53_sj_data_read[,"missed_exon_small_insh_count"]))/
#   as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);
# 
# g_53_sj_data_read[,"corrected missed_anno_small_exon_insh_count"]<-
#   as.numeric(as.character(g_53_sj_data_read[,"missed_anno_small_exon_insh_count"]))/
#   as.numeric(as.character(g_53_sj_data_read$read_count))*(10^6);



setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
library(ggplot2);
library(bedr)
options("scipen"=100, "digits"=4)

#       for i in `ls /Users/mengli/Documents/projects/abs/data/star_abs5_motif/*.seq`;do `perl score5.pl $i >> $i.score`;done;
#       for i in `ls /Users/mengli/Documents/projects/abs/data/star_abs3_motif/*.seq`;do `perl score3.pl $i >> $i.score`;done;

#       find . -type f -name "*SJ.out.tab" -exec cp {} ../star_re \;
#       find . -type f -name "*final.out" -exec cp {} ../star_log \;

#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/data/polII_mutant_data/star/\*SJ.out.tab \
#   /Users/mengli/Documents/projects/abs/data/poll_star/

#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/data/star_re/\* /Users/mengli/Documents/projects/abs/data/star/
#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/data/star_log/\* /Users/mengli/Documents/projects/abs/data/star_log/


#   find . -type d -name "tmp" -exec rm -rf {} \;


source("mouse/ana_star_full_mm_cgr8.R",echo = TRUE);
source("mouse/ana_star_full_mm_n2a.R",echo = TRUE);

g_53_sj_data<-rbind(g_53_sj_data_n2a,g_53_sj_data_cgr8);

write.table(g_53_sj_data,file="mouse/g_53_sj_all_full_mm.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);

gene_read_fr_fr<-read.table("mouse/rbp_read_count_map.tsv",header = TRUE,as.is = TRUE,sep="\t")


g_53_sj_data_read<-inner_join(g_53_sj_data,gene_read_fr_fr,by=c("rbp"="g") );

g_53_sj_data_read[,"corrected cryptic 3'ss"]<-as.numeric(as.character(g_53_sj_data_read[,"cryptic 3'ss"]))/
  as.numeric(as.character(g_53_sj_data_read$read_count) )*(10^6);

g_53_sj_data_read[,"corrected GC type 5'ss"]<-as.numeric(as.character(g_53_sj_data_read[,"GC type 5'ss"]))/
  as.numeric(as.character(g_53_sj_data_read$read_count) )*(10^6);

g_53_sj_data_read[,"is_N2A"]<-str_detect(g_53_sj_data_read$rbp,"N2A");

#g_53_sj_data_read$X3_X5_ratio<-g_53_sj_data_read$X3/g_53_sj_data_read$X5;

#colnames(g_53_sj_data_read)<-c("rbp","cryptic 5'ss","cryptic 3'ss",
#                               "cryptic whole junction","U12 cryptic 5'ss","U12 cryptic 3'ss","CG type 5'ss",
#                               "read_count","cryptic 3'ss correct read count","cryptic 3'ss/cryptic 5'ss");

write.table(g_53_sj_data_read,file="mouse/g_53_sj_all_read_count_mouse.tsv",sep="\t",
            quote = FALSE,col.names = TRUE,row.names = FALSE);


source("mouse/mouse_3ss_5ss_idr.R",echo = TRUE);




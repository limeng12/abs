# count the junction in each intron contain sno or not
library(stringr)
library(ggplot2)
library(GenomicFeatures)
setwd("/Users/mengli/Documents/projects/abs");

source("code/multiplot.R")
#intron_coor<-read.table(file="anno/intron_coor_gencode.tsv",sep="\t", as.is = TRUE,header = TRUE);

txdb <- makeTxDbFromGFF("anno/hg19_refgene.gtf");

intron_coor_range <- intronicParts(txdb);

#intron_coor_range<-with(intron_coor,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand) );



junc_coor<-read.table("data/star_target_only_jc/ENCSR624OUI_AQR_poly_no_ctl.bed",sep="\t",header  = TRUE,as.is = TRUE);
colnames(junc_coor)<-c("chr","start","end");

junc_coor_range<-with(junc_coor,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand="*") );



intron_overlap_number<-countOverlaps(intron_coor_range,junc_coor_range);

intron_within_number<-countOverlaps(intron_coor_range,junc_coor_range, type="within");

intron_overlap_jc_number<-intron_overlap_number-intron_within_number;



ref_gene<-read.table("anno/hg19_refgene.tsv",header = TRUE,sep = "\t",as.is = TRUE);

# SNORD, SNORA, SCARNA
ref_gene_sno<-ref_gene[str_detect(ref_gene$name2,"^SNORD") | 
                         str_detect(ref_gene$name2,"^SNORA") | 
                         str_detect(ref_gene$name2,"^SCARNA"),];



ref_gene_sno_gr<-with(ref_gene_sno ,GRanges(seqnames = chrom,IRanges(start=txStart,end=txEnd),strand=strand) );

intron_overlap_sno_number<-countOverlaps(intron_coor_range,ref_gene_sno_gr);


hk_genes<-read.table("anno/housekeeping_genes.tsv",header = FALSE,as.is = TRUE)

ref_gene<-read.table("anno/hg19_refgene.tsv",header = TRUE,sep = "\t",as.is = TRUE);

ref_gene_sel_hk<-ref_gene[ref_gene$name2 %in% hk_genes[,1],];


ref_gene_sel_hk_gr<-with(ref_gene_sel_hk ,GRanges(seqnames = chrom,IRanges(start=txStart,end=txEnd),strand=strand) );

if_hk_count<-countOverlaps(intron_coor_range,ref_gene_sel_hk_gr);

intron_len<-end(intron_coor_range)-start(intron_coor_range);

intron_len_cut<-(intron_len>3000);

jc_count_fr<-data.frame(jc_count=intron_overlap_jc_number,
                        sno_count=intron_overlap_sno_number,
                        intron_len=intron_len,
                        intron_len_cut=intron_len_cut,
                        if_hk=if_hk_count>0);


jc_count_fr_filter<-jc_count_fr[jc_count_fr$jc_count>0,];


jc_count_fr_filter$has_sno<-jc_count_fr_filter$sno_count>0;



p1<-ggplot(jc_count_fr_filter)+geom_boxplot(aes(x=has_sno,
                                              y=jc_count),fill="#E69F00" )+scale_y_log10()+
  theme_minimal()+theme(text = element_text(size=13))+xlab("")+ylab("")
  #ggtitle("Compared with all genes")+
  #xlab("if contain snoRNA")+
  #ylab("cryptic junction count")+ theme_minimal()+theme(text = element_text(size=13))


jc_count_fr_filter2<-jc_count_fr_filter[jc_count_fr_filter$sno_count>0 | jc_count_fr_filter$if_hk>0,];

p2<-ggplot(jc_count_fr_filter2)+geom_boxplot(aes(x=has_sno,
                                                y=jc_count),fill="#E69F00" )+scale_y_log10()+
  theme_minimal()+theme(text = element_text(size=13))+xlab("")+ylab("")
  #ggtitle("Compared with housekeeping genes")+
  #xlab("if contain snoRNA")+
  #ylab("cryptic junction count")+ 



p3<-ggplot(jc_count_fr_filter2[jc_count_fr_filter2$has_sno==TRUE,])+geom_boxplot(aes(x=intron_len_cut,
                                                 y=jc_count) )+scale_y_log10()+
  ggtitle("cryptic ss for long and short intron, both contain snoRNA")+
  xlab("intron length > 1000")+
  ylab("cryptic junction count");

jc_count_fr_filter2_sno_long<-jc_count_fr_filter2[ (jc_count_fr_filter2$has_sno==TRUE) & (jc_count_fr_filter2$intron_len_cut==TRUE),]


p4<-ggplot(jc_count_fr_filter2[jc_count_fr_filter2$has_sno==FALSE,])+geom_boxplot(aes(x=intron_len_cut,
                                                                                     y=jc_count) )+scale_y_log10()+
  ggtitle("cryptic ss for long and short intron, both not contain snoRNA")+
  xlab("intron length > 1000")+
  ylab("cryptic junction count");



p5<-ggplot(jc_count_fr_filter2[jc_count_fr_filter2$has_sno==TRUE,])+geom_point(aes(x=intron_len,
                                                                                     y=jc_count) )+scale_y_log10()+
  ggtitle("cryptic ss for different intron length contain snoRNA")+
  xlab("intron length")+
  ylab("cryptic junction count");

#jc_count_fr_filter2_sno_long<-jc_count_fr_filter2[ (jc_count_fr_filter2$has_sno==TRUE) & (jc_count_fr_filter2$intron_len_cut==TRUE),]


p6<-ggplot(jc_count_fr_filter2[jc_count_fr_filter2$has_sno==FALSE,])+geom_point(aes(x=intron_len,
                                                                                      y=jc_count) )+scale_y_log10()+
  ggtitle("cryptic ss for different intron length contain no snoRNA")+
  xlab("intron length ")+
  ylab("cryptic junction count");


pdf("result/jc_count_and_sno.pdf",width = 5,height = 2.5);

multiplot(p1,p2,cols=2);
dev.off();

#pdf("result/jc_count_and_sno.pdf",width = 16,height = 8);

#pdf("result/jc_count_and_sno.pdf",width = 16,height = 8);
#print(p1);


#print(p2);

#multiplot(p1,p3,p2,p4,p5,p6,cols=1);
#multiplot(p1,p3,p2,p4,cols=2);

#dev.off();


#intron_coor_junc_count<-cbind(intron_coor,intron_overlap_number);

#intron_junc_count_normal<-intron_overlap_number/(intron_coor$end-intron_coor$start);



#intron_coor_junc_count<-cbind(intron_coor,intron_overlap_number,intron_junc_count_normal);


#intron_coor_junc_count<-intron_coor_junc_count[
#  order(intron_coor_junc_count$intron_junc_count_normal,decreasing=TRUE),];

#write.table(intron_coor_junc_count,file="result/intron_coor_junc_count.tsv",
#            sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE);



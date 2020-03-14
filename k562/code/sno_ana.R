##count bs sites number in sno or not

library(GenomicFeatures);
library(stringr)
library(plyr)
library(dplyr)
library(GenomicRanges)

setwd("/Users/mengli/Documents/projects/abs");
source("code/multiplot.R")

txdb <- makeTxDbFromGFF("anno/hg19_refgene.gtf");

all.introns <- intronicParts(txdb);


elementMetadata(all.introns)[,"tx_name_one"]<-
  sapply(str_split(sapply(elementMetadata(all.introns)[,"tx_name"],"[",1),"\\."),"[",1);


intron_coor<-as.data.frame(cbind(as.character(seqnames(all.introns) ),start(all.introns),
                                 end(all.introns),as.character(strand(all.introns)) ) ,
                           stringsAsFactors=FALSE);

colnames(intron_coor)<-c("chr","start","end","strand");

intron_coor$start<-as.numeric(intron_coor$start);
intron_coor$end<-as.numeric(intron_coor$end);


intron_coor_gr<-with(intron_coor ,GRanges(seqnames = chr,IRanges(start=start,end=end),strand=strand) );


ref_gene<-read.table("anno/hg19_refgene.tsv",header = TRUE,sep = "\t",as.is = TRUE);

# SNORD, SNORA, SCARNA
ref_gene_sel<-ref_gene[str_detect(ref_gene$name2,"^SNORD") | 
                         str_detect(ref_gene$name2,"^SNORA") | 
                         str_detect(ref_gene$name2,"^SCARNA"),];

#ref_gene_sel<-ref_gene[,];


ref_gene_sel_gr<-with(ref_gene_sel ,GRanges(seqnames = chrom,IRanges(start=txStart,end=txEnd),strand=strand) );

if_sn<-countOverlaps(intron_coor_gr,ref_gene_sel_gr);

intron_coor<-cbind(intron_coor,if_sn);

elementMetadata(intron_coor_gr)[[ "if_sn" ]] <- if_sn


bp<-read.table("anno/FileS2_BS_site.bed",header = FALSE,as.is = FALSE,sep="\t");
colnames(bp)<-c("chr","start","end","V4","V5","strand")

bp_gr<-with(bp ,GRanges(seqnames = chr,IRanges(start=start,end=end),strand=strand) );


count_of_bp<-countOverlaps(intron_coor_gr,bp_gr);


elementMetadata(intron_coor_gr)[[ "count_of_bp" ]] <- count_of_bp


hk_genes<-read.table("anno/housekeeping_genes.tsv",header = FALSE,as.is = TRUE)

ref_gene<-read.table("anno/hg19_refgene.tsv",header = TRUE,sep = "\t",as.is = TRUE);

ref_gene_sel_hk<-ref_gene[ref_gene$name2 %in% hk_genes[,1],];


ref_gene_sel_hk_gr<-with(ref_gene_sel_hk ,GRanges(seqnames = chrom,IRanges(start=txStart,end=txEnd),strand=strand) );

if_hk<-countOverlaps(intron_coor_gr,ref_gene_sel_hk_gr);

elementMetadata(intron_coor_gr)[[ "is_hk" ]] <- (if_hk>0)


intron_coor_gr_fr<-as.data.frame(intron_coor_gr);


write.table(intron_coor_gr_fr,file="result/intron_coor_gr_fr.tsv",
            sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)


intron_coor_gr_fr$`has sno`<-as.factor(intron_coor_gr_fr$if_sn>0)

library(plyr)
intron_coor_gr_fr$`has sno`<-revalue(intron_coor_gr_fr$`has sno`, c("FALSE"="don't have sno", "TRUE"="have sno"))


#levels(as.factor(intron_coor_gr_fr$`has sno`))<-c("don't have sno","have sno")


intron_coor_gr_fr_sel<-intron_coor_gr_fr[intron_coor_gr_fr$`has sno`=="have sno" | intron_coor_gr_fr$is_hk,]


#t.test(intron_coor_gr_fr[intron_coor_gr_fr$if_sn>0,"count_of_bp"],
#       intron_coor_gr_fr[intron_coor_gr_fr$if_sn==0,"count_of_bp"])



#boxplot(intron_coor_gr_fr[intron_coor_gr_fr$if_sn>0,"count_of_bp"],
#        intron_coor_gr_fr[intron_coor_gr_fr$if_sn==0,"count_of_bp"], log="y")

library(ggplot2)
# Basic box plot
p1 <- ggplot(intron_coor_gr_fr, aes(x=`has sno`, y=count_of_bp)) + scale_y_log10()+
  geom_boxplot()+xlab("compared with all genes")+theme_minimal()+ theme(text = element_text(size=20) )+
  ylab("# of branch points")+ggtitle("# of branch points of normal intron and snoRNA contained intron")+theme_minimal()

p2 <- ggplot(intron_coor_gr_fr_sel, aes(x=`has sno`, y=count_of_bp)) + scale_y_log10()+
  geom_boxplot()+xlab("compared with only housekeeping genes")+theme_minimal()+ theme(text = element_text(size=20) )+
  ylab("# of branch points")+theme_minimal()



pdf("result/sno_compared_all.pdf", width=12,height=6)
#print(p1)

#print(p2)
multiplot(p1,p2,cols = 2)
dev.off();




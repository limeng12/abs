## maxEntScore for each cryptic 5ss,3ss and intron
setwd("/Users/mengli/Documents/projects/abs");

library(BSgenome.Hsapiens.UCSC.hg19);
library(GenomicFeatures)
library(stringr)
library(plyr)
library(dplyr)
library(RWebLogo)

hg19_genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19");
setwd("/Users/mengli/Documents/projects/abs");



txdb <- makeTxDbFromGFF("anno/hg19_refgene.gtf");

all.introns <- intronicParts(txdb);

all.introns_cryptic_intron<-all.introns[!duplicated(all.introns),];

intron_coor<-data.frame(chr=as.character(seqnames(all.introns_cryptic_intron) ),
                        start=start(all.introns_cryptic_intron),
                        end=end(all.introns_cryptic_intron),
                        strand=as.character(strand(all.introns_cryptic_intron)) ,
                        stringsAsFactors=FALSE);

intron_coor_unique<-unique(intron_coor);

intron_coor_positive<-intron_coor_unique[intron_coor_unique$strand=="+",];


t_5_seqs_sel<-t_5_seqs<-getSeq(hg19_genome,intron_coor_positive[,"chr"],
                 intron_coor_positive[,"start"]-3,
                 intron_coor_positive[,"start"]+5,
                 strand=intron_coor_positive[,"strand"],as.character=TRUE);

#t_5_seqs_sel<-t_5_seqs[str_sub(t_5_seqs,4,5)=="GT"];

cat(as.character(t_5_seqs_sel),file=paste0("anno/intron5.seq"), sep="\n");

weblogo(t_5_seqs_sel,file.out=paste0("result/intron5_motif.pdf") ,open = FALSE );




t_3_seqs_sel<-t_3_seqs<-getSeq(hg19_genome,intron_coor_positive[,"chr"],
                 intron_coor_positive[,"end"]-19, 
                 intron_coor_positive[,"end"]+3,
                 strand=intron_coor_positive[,"strand"],as.character=TRUE);

#t_3_seqs_sel<-t_3_seqs[str_sub(t_3_seqs,19,20)=="AG"];

cat(as.character(t_3_seqs_sel),file=paste0("anno/intron3.seq"), sep="\n");

weblogo(t_3_seqs_sel,file.out=paste0("result/intron3.pdf") ,open = FALSE );




setwd("/Users/mengli/Documents/software/maxentscan");

system("perl score5.pl /Users/mengli/Documents/projects/abs/anno/intron5.seq > /Users/mengli/Documents/projects/abs/anno/intron5.seq.score")

system("perl score3.pl /Users/mengli/Documents/projects/abs/anno/intron3.seq > /Users/mengli/Documents/projects/abs/anno/intron3.seq.score")


system("perl score5.pl /Users/mengli/Documents/projects/abs/data/star_abs5_motif/ENCSR624OUI_AQR_poly.seq > /Users/mengli/Documents/projects/abs/data/star_abs5_motif/ENCSR624OUI_AQR_poly.seq.score")

system("perl score3.pl /Users/mengli/Documents/projects/abs/data/star_abs3_motif/ENCSR624OUI_AQR_poly.seq > /Users/mengli/Documents/projects/abs/data/star_abs3_motif/ENCSR624OUI_AQR_poly.seq.score")


setwd("/Users/mengli/Documents/projects/abs");



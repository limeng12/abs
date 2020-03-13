## score and motif of all AG in introns

library(GenomicFeatures)
library(stringr)
library(plyr)
library(dplyr)
library(RWebLogo)

setwd("/Users/mengli/Documents/projects/abs");


txdb <- makeTxDbFromGFF("anno/hg19_refgene.gtf");

all.introns <- intronicParts(txdb);


elementMetadata(all.introns)[,"tx_name_one"]<-
  sapply(str_split(sapply(elementMetadata(all.introns)[,"tx_name"],"[",1),"\\."),"[",1);


intron_coor<-as.data.frame(cbind(as.character(seqnames(all.introns) ),start(all.introns),
                                 end(all.introns),as.character(strand(all.introns)) ) ,
                           stringsAsFactors=FALSE);

colnames(intron_coor)<-c("chr","start","end","strand")

###only consider positive strand
intron_coor<-intron_coor[intron_coor$strand=="+",]
  
intron_coor$start<-as.numeric(intron_coor$start);
intron_coor$end<-as.numeric(intron_coor$end);

intron_coor_gr<-with(intron_coor ,
                     GRanges(seqnames = chr,IRanges(start=start,end=end),strand=strand) );


#intron_coor_gr_reduce<-reduce(intron_coor_gr);
###########delete super big introns that contain two or more small introns. 
#count_over<-countOverlaps(intron_coor_gr,intron_coor_gr,type="within");
#intron_coor_gr_sel<-intron_coor_gr[count_over<=1,];


library("BSgenome");
library(BSgenome.Hsapiens.UCSC.hg19);
genome <- BSgenome.Hsapiens.UCSC.hg19;

#getSeq(hg38_genome,data_filter[,"Region"],data_filter[,"Position"]-49,data_filter[,"Position"]+50);

intron_seq_all<-getSeq(genome,intron_coor_gr);

seq_list_all<-list();

seq_list_all<-sapply(intron_seq_all,as.character);


seq_list_all_v<-as.character(seq_list_all);

AG_pos<-str_locate_all(seq_list_all_v,"AG");

intron_start<-start(intron_coor_gr_sel);

intron_end<-end(intron_coor_gr_sel);

intron_seqnames<-seqnames(intron_coor_gr_sel);

all_ag_seq<-c();
for(i in 1:length(AG_pos) ){
  one_intron_ag<-AG_pos[[i]];
  
  if(nrow(one_intron_ag)<2){
    next;
  }
  
  for(j in 1:(nrow(one_intron_ag)-1) ){
    if(length(all_ag_seq)>10000){
      break;
    }
    
    #print(paste0("start pos: ",intron_start[i]+one_intron_ag[j,"start"]," chr: ",intron_seqnames[i]) );
    
    aaaa<-as.character(  getSeq(genome,intron_seqnames[i],intron_start[i]+one_intron_ag[j,"start"]-19,
          intron_start[i]+one_intron_ag[j,"start"]+3) );
    
    all_ag_seq<-c(all_ag_seq,aaaa);
    
  }
  
}


weblogo(all_ag_seq,file.out=paste0("result/anno_all_cryptic_ag_score.pdf") ,open = FALSE )


writeLines(all_ag_seq,con = file("result/anno_all_cryptic_ag_score.seq") );

setwd("/Users/mengli/Documents/software/maxentscan");

system(" perl score3.pl /Users/mengli/Documents/projects/abs/result/anno_all_cryptic_ag_score.seq > \
       /Users/mengli/Documents/projects/abs/result/anno_all_cryptic_ag_score.score")


setwd("/Users/mengli/Documents/projects/abs");





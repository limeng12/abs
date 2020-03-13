library(ggplot2)
library(stringr)

setwd("/Users/mengli/Documents/projects/abs");

####cryptic site score

data3_mat<-matrix(nrow=0,ncol=2);


files_all<-list.files("/Users/mengli/Documents/projects/abs/data/star_abs5_motif","*.score");

gene_ids<-sapply(str_split(files_all,"\\.seq"),"[",1);

#gene_ids<-"ENCSR624OUI_AQR_poly";

for( i in gene_ids){
  file_path<-paste0("/Users/mengli/Documents/projects/abs/data/star_abs5_motif/",i,".seq.score");
  
  value_fr<-read.table(file_path,sep="\t",header = TRUE,as.is = TRUE);
  names<-rep("used_cryptic",nrow(value_fr));
  
  data3_mat<-rbind(data3_mat,cbind(names,value_fr[,2] ) );
  
}

data3_mat<-data3_mat[sample(1:nrow(data3_mat),100000),];


###intron score
value_fr<-read.table("anno/intron5.seq.score",sep="\t",header = TRUE,as.is = TRUE);
names<-rep("intron",nrow(value_fr));
median(value_fr[,2]);#8.48

data3_mat<-rbind(data3_mat,cbind(names,value_fr[,2] ) );





value_fr2<-read.table("result/anno_all_cryptic_gt_score.score",sep="\t",header = TRUE,as.is = TRUE);
names<-rep("all_intron_cryptic",nrow(value_fr2));

data3_mat<-rbind(data3_mat,cbind(names,value_fr2[,2] ) );

median(value_fr2[,2])

data3_mat_fr<-data.frame(name=as.factor(data3_mat[,1]),
                         value=as.numeric(data3_mat[,2]) );


p <- ggplot(data3_mat_fr, aes(x=name, y=value)) + 
  geom_boxplot(fill="#E69F00")+theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("maxEntScore")+theme_minimal()+
  scale_x_discrete(limits=c("all_intron_cryptic","used_cryptic","intron"),
                   labels=c("All cryptic 5' sites"," Used cryptic 5' sites","Canonical intron 5' splice site"))+
  xlab("")#+ggtitle("maxEntScore of different 5' splice sites")



#source("code/multiplot.R");
pdf("result/cryptic_ss_vs_real_ss_VS_all_crytic_5ss.pdf",width=7.5,height=5)
print(p)
dev.off()



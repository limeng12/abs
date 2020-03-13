# cluster based on RBP shared junctions, may induce regulation or a complex

setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);
library(cluster);
library(pheatmap)
library(reshape2)



gene_ids<-(unique(readLines("mouse/gene_id_mm.txt")) );

gene_ids_noctl<-gene_ids[!str_detect(gene_ids,"siNT|mock|untreated")];

jc_list<-list();

g_no_input<-c();

all_jc<-c()

for(g in gene_ids_noctl){
  
  if(!file.exists(paste0("mouse/star_target_only_jc/",g,"_no_ctl.bed") )){
    g_no_input<-c(g_no_input,g)
    next;
  }
  
  if(file.size(paste0("mouse/star_target_only_jc/",g,"_no_ctl.bed") )==0){
    g_no_input<-c(g_no_input,g)
    
    next;
  }
  
  
  t_sj<-read.table(paste0("mouse/star_target_only_jc/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE)[,1:7];
  
  if(nrow(t_sj)<10){
    g_no_input<-c(g_no_input,g)
    
    next;
  }
  
  
  colnames(t_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n");
  
  
  all_jc<-c(all_jc,str_c(t_sj[,"chr"],":",t_sj[,"X5_pos"],"-",t_sj[,"X3_pos"]) );
  
  jc_list[[g]]<-str_c(t_sj$chr,":",t_sj$X5_pos,"-",t_sj$X3_pos);
  
}

gene_ids_noctl<-setdiff(gene_ids_noctl,g_no_input);

diss_mat<-matrix(nrow=length(gene_ids_noctl),ncol=length(gene_ids_noctl) );

rownames(diss_mat)<-gene_ids_noctl;
colnames(diss_mat)<-gene_ids_noctl;

#jc_total_tbl<-read.table("result/all_junc_nodup.tsv",sep = "\t", header = TRUE,as.is = TRUE);

#jc_total<-nrow(jc_total_tbl);

jc_total<-length(unique(all_jc));


cat("RBP1\tRBP2\tvalue\tA and B\tA\tB\tA minus B\t B minus A\tnot in AB\n",
    file="mouse/jc_overlap_p_values_melt_mm.tsv");

for(g1 in gene_ids_noctl){
  
  for(g2 in gene_ids_noctl){
    
    a<-length(jc_list[[g1]]);
    
    b<-length(jc_list[[g2]]);
    
    
    aa<-length(intersect(jc_list[[g1]],jc_list[[g2]]) )
    
    ab<-length(setdiff(jc_list[[g1]],jc_list[[g2]]));
    
    ba<-length(setdiff(jc_list[[g2]],jc_list[[g1]]));
    
    aabb<-jc_total-length(union(jc_list[[g2]],jc_list[[g1]]));
    
    f_mat<-matrix(c(aa,
                    ab,
                    ba,
                    aabb
                    ),nrow=2,byrow = TRUE);
    
    #diss_mat[g1,g2]<- (fisher.test(f_mat))$p.value;
    
   # exp_aa<-(ab/aabb*ba);
    exp_aa<-(a/jc_total*b)
    
    if(g1!=g2){
      diss_mat[g1,g2]<-ppois(aa,exp_aa,lower.tail = FALSE,log.p = TRUE);
    }else{
      diss_mat[g1,g2]<-NA;
    }
    
    
    if(g1>g2){
      
      cat(paste0(g1,"\t",g2,"\t",diss_mat[g1,g2],"\t",aa,"\t",a,"\t",b,"\t",ab,"\t",ba,"\t",aabb,"\n"),
          sep="",file="mouse/jc_overlap_p_values_melt_mm.tsv",append = TRUE);
      
    }
    
    #diss_mat[g1,g2]<-length(intersect(jc_list[[g1]],jc_list[[g2]]) );
     
  }
}

rownames(diss_mat)<-gene_ids_noctl;
colnames(diss_mat)<-gene_ids_noctl;



write.table(diss_mat,file = "mouse/jc_overlap_p_values_matrix_mm.tsv",sep = "\t",
            col.names = TRUE,row.names = TRUE,quote = FALSE);

library(reshape2)


diss_mat<- log(-1*diss_mat+1);


pdf("mouse/jc_ovlap_clu_mm.pdf",width = 100,height = 100);

pheatmap(diss_mat,cluster_rows=TRUE,cluster_cols=TRUE);
#pheatmap(diss_mat,cluster_rows=TRUE,cluster_cols=TRUE,breaks=unique(c(0,0.92,seq(0.92,1,by=0.001) )),
#         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(81)  );

dev.off();

# color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

#LIN28B NUSAP1 
#132    208 




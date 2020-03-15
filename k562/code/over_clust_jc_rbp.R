setwd("/Users/mengli/Documents/projects/abs");
library(plyr);
library(stringr);
library(dplyr);
library(cluster);
library(pheatmap)
library(reshape2)

# gene_ids<-readLines("samples/gene_id.txt");
gene_ids_sicr<-(unique(readLines("samples/gene_id_sicr.txt") ) );
gene_ids_bren<-(unique(readLines("samples/gene_id.txt") ) );

gene_ids<-c(gene_ids_bren,gene_ids_sicr)#,gene_ids_sicr);
gene_ids<-gene_ids[!str_detect(gene_ids,"_inte")]

gene_ids_noctl<-gene_ids[!str_detect(gene_ids,"CTL_")];

all_jc<-c()

jc_list<-list();

for(g in gene_ids_noctl){
  
  t_sj<-read.table(paste0("data/star_target_only_jc/",g,"_no_ctl.bed"),sep = "\t",header = FALSE, as.is = TRUE)[,1:7];
  
  colnames(t_sj)<-c("chr","X5_pos","X3_pos","strand","type","X5_n","X3_n");
  
  jc_list[[g]]<-str_c(t_sj$chr,":",t_sj$X5_pos,"-",t_sj$X3_pos);
  
  
  all_jc<-c(all_jc,str_c(t_sj[,"chr"],":",t_sj[,"X5_pos"],"-",t_sj[,"X3_pos"]) );
  
}


diss_mat<-matrix(nrow=length(gene_ids_noctl),ncol=length(gene_ids_noctl) );

rownames(diss_mat)<-gene_ids_noctl;
colnames(diss_mat)<-gene_ids_noctl;

#jc_total_tbl<-read.table("result/all_junc_nodup.tsv",sep = "\t", header = TRUE,as.is = TRUE);

#jc_total<-nrow(jc_total_tbl);

jc_total<-length(unique(all_jc));

cat("RBP1\tRBP2\tvalue\tA and B\tA\tB\tA minus B\t B minus A\tnot in AB\n",
    file="result/jc_overlap_p_values_melt.tsv");

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
    #exp_aa<-(ab/aabb*ba)
    exp_aa<-(a/jc_total*b);
    
    
    if(g1!=g2){
      diss_mat[g1,g2]<-ppois(aa,exp_aa,lower.tail = FALSE,log.p = TRUE);
    }else{
      diss_mat[g1,g2]<-NA;
    }
    
    if(g1>g2){
      
      cat(paste0(g1,"\t",g2,"\t",diss_mat[g1,g2],"\t",aa,"\t",a,"\t",b,"\t",ab,"\t",ba,"\t",aabb,"\n"),
          sep="",file="result/jc_overlap_p_values_melt.tsv",append = TRUE);
      
    }
    
    #diss_mat[g1,g2]<-length(intersect(jc_list[[g1]],jc_list[[g2]]) );
     
  }
}

rownames(diss_mat)<-gene_ids_noctl;
colnames(diss_mat)<-gene_ids_noctl;


#write.table(diss_mat,file = "result/jc_overlap_p_values_matrix.tsv",sep = "\t",
#            col.names = TRUE,row.names = TRUE,quote = FALSE);


library(reshape2)

diss_mat<- log(-1*diss_mat+1);

pdf("result/jc_ovlap_clu.pdf",width = 50,height = 50);

#pheatmap(diss_mat,cluster_rows=TRUE,cluster_cols=TRUE);
pheatmap(diss_mat,cluster_rows=TRUE,cluster_cols=TRUE);

#pheatmap(diss_mat,cluster_rows=TRUE,cluster_cols=TRUE,breaks=unique(c(0,0.95,seq(0.95,1,by=0.001) )),
#         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(52)  );

dev.off();



#LIN28B NUSAP1 
#132    208 

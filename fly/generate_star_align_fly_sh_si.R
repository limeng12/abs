library(stringr);


#   scp limeng@10.10.118.191:/picb/rnasys2/limeng/fly/star_fly_sorted/ENCSR566KSW_Rox8/\*.bam /Volumes/mengli/abs/fly/bam/


# rm(list=ls())
# clean mapsplice bams
#  setsid find . -type f -name "*.bam" -exec rm -f {} \;
#  setsid find . -type d -name "tmp" -exec rm -rf {} \;
#  setsid find . -type d -name "tmp" -exec rm -rf {} +
#  find . -type f -name "*.bak" -exec rm -f {} \;


#  CIRCexplorer2 parse -t STAR -b aaa ZRANB2_1Chimeric.out.junction

#  sh circular_run.sh
#  sh circular_rename.sh

setwd("/Users/mengli/Documents/projects/abs");

#data_map<-read.table("fly/samples/fly_siRNA_seq_meta.tsv",header = TRUE ,sep = "\t", as.is = TRUE,comment.char="$");
data_map1<-read.table("fly/samples/fly_embryo_siRNA_all_meta.tsv",header = TRUE ,sep = "\t", as.is = TRUE,comment.char="$");
data_map2<-read.table("fly/samples/fly_shRNA_susan_metadata.tsv",header = TRUE ,sep = "\t", as.is = TRUE,comment.char="$");
data_map<-rbind(data_map1,data_map2);

data_map<-data_map2

#data_map_filter<-data_map[nchar(data_map$r1_q1)>5,];

#data_map_filter<-data_map[(data_map[,"Library.made.from"]=="polyadenylated mRNA") |
#                            (data_map[,"Experiment.target"]=="Non-specific target control-human"),];


data_map_fq<-data_map[data_map$File.format=="fastq" &
                        (data_map$Biosample.treatments=="" | is.na(data_map$Biosample.treatments)) ,];
#data_map_filter<-data_map_filter[!duplicated(data_map_filter$Target.gene),];

#fastq_ids<-c(data_map_filter$r1_q1,data_map_filter$r1_q2,data_map_filter$r2_q1,data_map_filter$r2_q2)

fastq_ids<-data_map_fq$File.download.URL;
fastq_ids<-str_c("",fastq_ids);

#unlink("fly/download_fastq.sh");

#cat(c("#!/bin/bash",fastq_ids),file="fly/download_fastq.sh",sep = "\n",append = TRUE);
#cat(c("#!/bin/bash",fastq_ids),file="fly/download_fastq.sh",sep = "\n",append = TRUE);

data_map_filter<-data_map_fq;

data_map_filter$Experiment.target<-str_sub(data_map_filter$Experiment.target,1,-15)


#data_map_filter_onepair<-data_map_filter[data_map_filter$Paired.end==1,];


unlink("fly/code/run_sh/align_fly_all.sh");
#unlink("fly/code/run_sh/align_fly_sorted.sh");

#unlink("fly/code/run_sh/cufflink_fly.sh");
#unlink("fly/code/run_sh/circular_run_fly.sh");#circular_rename.sh
#unlink("fly/code/run_sh/circular_rename_fly.sh");
unlink("fly/samples/gene_id_fly_all.txt");
#unlink("fly/code/run_sh/midec_fly.sh");


cat("STAR --genomeLoad Remove --genomeDir /picb/rnasys2/limeng/anno/dm6/ENSEMBL/star_index\n",
    file="fly/code/run_sh/align_fly_all.sh",append = TRUE);

cat("STAR --genomeLoad LoadAndExit --genomeDir /picb/rnasys2/limeng/anno/dm6/ENSEMBL/star_index\n\n",
    file="fly/code/run_sh/align_fly_all.sh",append = TRUE);

cat("rm star_wrong.txt\n",file = "fly/code/run_sh/align_fly_all.sh",append = TRUE);

t<-1;

control_index<-0;

#data_map_filter_onepair<-
#  data_map_filter_onepair[with(data_map_filter_onepair,
#                               order(Experiment.accession,Biological.replicate.s.)),];

data_map_filter<-
  data_map_filter[!duplicated(data_map_filter[,c("Experiment.accession","Biological.replicate.s.","File.accession")]),];

#data_map_filter_onepair<- data_map_filter_onepair[seq(nrow(data_map_filter_onepair),1),];


for(i in 1:nrow(data_map_filter) ){
  gene<-data_map_filter$Experiment.target[i];
  
  
  exp_accession<-data_map_filter$Experiment.accession[i];
  
  file.accession<-data_map_filter$File.accession[i];
  
  assay<-str_sub(data_map_filter$Assay[i],1,5);
  
  
  rep_index<-data_map_filter$Biological.replicate.s.[i];
  #if(!(gene %in% splicing_rbps ) && !(str_detect(gene,"CTL_"))){
  #  next;
  #}
  
  if(str_detect(gene,"Control") | gene=="" ){
    if(rep_index==1){
      control_index<-control_index+1;
    }
    #gene<-paste0("CTL_",exp_accession,,"_",control_index);
    #gene<-paste0("CTL_",exp_accession,"_");
    gene<-paste0(exp_accession,"_CTL_",assay);
    
    
  }else{
    gene<-paste0(exp_accession,"_",gene,"_",assay);
    
  }
  
  
  #data_map_filter$Experiment.target[i]<-gene;
  
  cat(paste0(gene,"\n"),file="fly/samples/gene_id_fly_all.txt",append = TRUE);
  
  #if(rep_index%%2==0){
  #  cat(paste0(gene,"\n"),file="samples/gene_id.txt",append = TRUE);
  #}
  
  #bam_name<-paste0("star/",gene,"/",gene,"_",rep_index,"Aligned.sortedByCoord.out.bam\n");
  #bam_name2<-paste0("star/",gene,"/",gene,"_2Aligned.sortedByCoord.out.bam\n");
  #cat(bam_name,file="code/run_sh/bamfile_names",append = TRUE);
  
  
  r1_q1<-data_map_filter$File.accession[i];
  #r1_q2<-data_map_filter$Paired.with[i];
  #r2_q1<-data_map_filter[i,3];
  #r2_q2<-data_map_filter[i,4];  
  #if(r1_q1=="ENCFF043USC"){
  #  aaa<-1
  #}
  cat(paste0("echo ", gene," >> star_wrong.txt\n"),file="fly/code/run_sh/align_fly_all.sh",append = TRUE);
  
  #align<-paste0("STAR --runThreadN 18 --quantMode GeneCounts --genomeLoad LoadAndKeep \\
  align<-paste0("STAR --runThreadN 18 --twopassMode Basic --twopass1readsN -1 --quantMode GeneCounts --genomeLoad NoSharedMemory \\
                --genomeDir /picb/rnasys2/limeng/anno/dm6/ENSEMBL/star_index \\
                --readFilesIn ../fly/fly_fastq/",r1_q1,".fastq.gz \\
                --chimJunctionOverhangMin 20 --chimSegmentMin 20 --chimOutType Junctions --chimSegmentReadGapMax 3 \\
                --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \\
                --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999  --bamRemoveDuplicatesType UniqueIdentical \\
                --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --readFilesCommand zcat \\
                --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitBAMsortRAM 30000000000 --outBAMsortingThreadN 16 \\
                --outSAMattributes NH HI NM MD AS nM jM jI XS --outFileNamePrefix \\
                ./star_fly/",gene,"/",gene,"_",file.accession,"_",rep_index," --outSAMtype BAM Unsorted 2>&1 | tee -a star_wrong.txt \n");
  #--chimJunctionOverhangMin 20 --chimSegmentMin 20 --chimOutType Junctions --chimSegmentReadGapMax 3 \\
  
  #if(rep_index==1){
  cat(paste0("mkdir ./star_fly/",gene,"\n"),file="fly/code/run_sh/align_fly_all.sh",append = TRUE);
  
  #}
  
  align_sort<-paste0("STAR --runThreadN 18 --twopassMode Basic --twopass1readsN -1 --quantMode GeneCounts --genomeLoad NoSharedMemory \\
                     --genomeDir /picb/rnasys2/limeng/anno/dm6/ENSEMBL/star_index \\
                     --readFilesIn ../fly/fly_fastq/",r1_q1,".fastq.gz \\
                     --chimJunctionOverhangMin 20 --chimSegmentMin 20 --chimOutType Junctions --chimSegmentReadGapMax 3 \\
                     --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \\
                     --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999  --bamRemoveDuplicatesType UniqueIdentical \\
                     --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --readFilesCommand zcat \\
                     --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitBAMsortRAM 30000000000 --outBAMsortingThreadN 16 \\
                     --outSAMattributes NH HI NM MD AS nM jM jI XS --outFileNamePrefix \\
                     ./star_fly_sorted/",gene,"/",gene,"_",file.accession,"_",rep_index," --outSAMtype BAM Unsorted 2>&1 | tee -a star_wrong.txt \n");
  
  
  #cat(paste0("mkdir ./star_fly_sorted/",gene,"\n"),file="fly/code/run_sh/align_fly_sorted.sh",append = TRUE);
  
  
  
  # cat(paste0("/picb/rnasys/program/install/java/jdk1.8.0_45/bin/java -jar midec.jar exon/",
  #            gene,"_",file.accession,"_",rep_index,"_exon star_fly/",gene,"/",gene,"_",file.accession,"_",rep_index,"Aligned.out.bam\n"),
  #     file="fly/code/run_sh/midec_fly.sh",append = TRUE); 
  
  
  #  cat(paste0("fast_circ.py parse -r /picb/rnasys2/limeng/anno/dm6/ENSEMBL/dm6_refgene_tbl_cir.txt -g /picb/rnasys2/limeng/anno/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa \\
  #             -t STAR -o circular_fly/",gene,"_",rep_index,"_circular star_fly/",gene,"/",gene,"_",file.accession,"_",rep_index,"Chimeric.out.junction\n"),
  #      file="fly/code/run_sh/circular_run_fly.sh",append = TRUE); 
  
  # cat(paste0("fast_circ.py parse -r /picb/rnasys2/limeng/anno/dm6/ENSEMBL/dm6_ensembl_tbl.txt -g /picb/rnasys2/limeng/anno/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa \\
  #            -t STAR -o circular_fly/",gene,"_",rep_index,"_circular star_fly/",gene,"/",gene,"_",file.accession,"_",rep_index,"Chimeric.out.junction\n"),
  #     file="fly/code/run_sh/circular_run_fly.sh",append = TRUE); 
  # 
  # cat(paste0("cp circular_fly/",gene,"_",rep_index,"_circular/circularRNA_known.txt circular_fly/",gene,"_",file.accession,"_",rep_index,"circularRNA_known.txt\n"),
  #     file="fly/code/run_sh/circular_rename_fly.sh",append = TRUE); 
  # cat(paste0("cp circular_fly/",gene,"_",rep_index,"_circular/back_spliced_junction.bed circular_fly/",gene,"_",file.accession,"_",rep_index,"back_spliced_junction.bed\n"),
  #     file="fly/code/run_sh/circular_rename_fly.sh",append = TRUE); 
  # 
  
  #cat(paste0("mkdir ./cufflink_fly/",gene,"_",file.accession,"_",rep_index,"\n"),file="code/run_sh/cufflink_fly.sh",append = TRUE);
  
  #link<-paste0("cufflinks -g /picb/rnasys/share/database/Homo_sapiens/GENCODE/hg19/gencode.v28lift37.annotation.gtf \\
  #             -F 0.0000001 -j 0.0000001 --min-frags-per-transfrag 1 -o cufflink_fly/",gene,"_",file.accession,"_",rep_index," -p 16 \\
  #             star_fly/",gene,"/",gene,"_",file.accession,"_",rep_index,"Aligned.sortedByCoord.out.bam 2>&1 | tee -a cufflink_fly_wrong.txt \n");
  
  #cat(link,file="code/run_sh/cufflink_fly.sh",append = TRUE);
  #cat(paste0("mkdir ./cuffmerge/",gene,"\n"),file="code/run_sh/cufflink_fly.sh",append = TRUE);
  
  
  cat(align,file="fly/code/run_sh/align_fly_all.sh",append = TRUE);
  
  #cat(align_sort,file="fly/code/run_sh/align_fly_sorted.sh",append = TRUE);
  
}


cat("STAR --genomeLoad Remove --genomeDir /picb/rnasys2/limeng/anno/dm6/ENSEMBL/star_index\n",
    file="fly/code/run_sh/align_fly_all.sh",append = TRUE);

#cat("normal",file="samples/gene_id.txt",append = TRUE);

system("scp /Users/mengli/Documents/projects/abs/fly/code/run_sh/align_fly_all.sh limeng@10.10.118.191:/picb/rnasys2/limeng/fly")


ids<-unique(readLines("fly/samples/gene_id_fly_all.txt") );
writeLines(ids,con=file("fly/samples/gene_id_fly_all.txt"))



#system("scp /Users/mengli/Documents/projects/abs/fly/code/run_sh/align_fly_sorted.sh limeng@10.10.118.191:/picb/rnasys2/limeng/fly")

#system("scp /Users/mengli/Documents/projects/abs/fly/code/run_sh/midec_fly.sh limeng@10.10.118.191:/picb/rnasys2/limeng/fly")

#system("scp /Users/mengli/Documents/projects/abs/fly/code/run_sh/circular_run_fly.sh limeng@10.10.118.191:/picb/rnasys2/limeng/fly")

#system("scp /Users/mengli/Documents/projects/abs/fly/code/run_sh/circular_rename_fly.sh limeng@10.10.118.191:/picb/rnasys2/limeng/fly")


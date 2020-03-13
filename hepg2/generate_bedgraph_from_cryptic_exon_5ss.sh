cd /Users/mengli/Documents/projects/abs/hepg2/data

#for i in `ls exon_mer_target_only_small/*_start.bed | xargs -n 1 basename`; \
#do `awk '{print $1"\t"$2"\t"$2+1"\t"$5}' exon_mer_target_only_small/$i > exon_mer_target_only_bedgraph/start_$i.bedgraph`;done;      
#for i in `ls exon_mer_target_only_small/*_end.bed | xargs -n 1 basename`; \
#do `awk '{print $1"\t"$3"\t"$3+1"\t"$5}' exon_mer_target_only_small/$i > exon_mer_target_only_bedgraph/end_$i.bedgraph`;done;   

#rm exon_mer_target_only_bedgraph/*
#  for i in `ls exon_mer_target_only_5ss/*.bed | xargs -n 1 basename`; \
#do `awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$5}' exon_mer_target_only_5ss/$i > exon_mer_target_only_bedgraph/$i.bedgraph`;done;      


rm exon_mer_target_only_bedgraph/*
  for i in `ls exon_mer_target_only_5ss/*.bed | xargs -n 1 basename`; \
do `awk '{print $1"\t"$2"\t"$3"\t"$6}' exon_mer_target_only_5ss/$i > exon_mer_target_only_bedgraph/$i.bedgraph`;done;      


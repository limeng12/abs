cd /Users/mengli/Documents/projects/abs/hepg2/data

for i in `ls star_target_only_jc_5ss/*.bed | xargs -n 1 basename`; \
do `awk '{print $1"\t"$2"\t"$2+1"\t"$6}' star_target_only_jc_5ss/$i > star_target_only_jc_bedgraph/5_$i.bedgraph`;done;  

for i in `ls star_target_only_jc_3ss/*.bed | xargs -n 1 basename`; \
do `awk '{print $1"\t"$2"\t"$2+1"\t"$6}' star_target_only_jc_3ss/$i > star_target_only_jc_bedgraph/3_$i.bedgraph`;done;   


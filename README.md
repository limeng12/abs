## Analysis of cryptic splice sites, cryptic juntions and cryptic exons.

Since the whole dataset is too big to upload into Github, this repository includes intermediate results and source code for this work.  

The result, source code and data for each organism/cell line are stored seperatelly in each folder. 

`K562` stores source code, motifs and cryptic splice sites related to K562 cell line in humans.

`HepG2` stores source code, motifs and cryptic splice sites related to HepG2 cell line in humans.

`mouse` stores source code, motifs and cryptic splice sites related to N2A and CGR8 cell lines in mice.

`fly` stores source code, motifs and cryptic splice sites related to S2 cell line in fruit flies.


### Source code (.R)
 `code` stores source codes.


### Cryptic splice sites and cryptic junctions (.bed)
`STAR` aligner will output `SJ` files, which is very convient for splice junction analysis. 

Read count support cryptic splice sites account 5/10,000 of total mapped reads in average, this make sensitivity is key thing to consider when analyzing cryptic splice sites. Thus I choose `STAR` with 2-pass mode to do alignment since it can detect highest number of cryptic splice sites. the sub-folders that starts with `star` are related to cryptic splice sites.

 `star_target_only_jc` stores the position of cryptic junctions.

 `star_target_only_jc_5ss` stores the position of cryptic 5' splice sites.

 `star_target_only_jc_3ss` stores the position of cryptic 3' splice sites.

### cryptic exons (.bed)

The `jI` tag in the `STAR` aligner, contain the junction sites, the regions between the junction sites are either introns or exons. I used this method to extract exon regions. 

 `exon_mer_target_only` stores the postion of cryptic exons.

 `exon_mer_target_only_5ss` stores the position of cryptic exons' 5' splice sites.

### Motifs of cryptic splie sites in different regions (.pdf)

I store cryptic splice sites' motifs of different regions in different folders. 

 `star_abs3_motif` stores the motifs of cryptic splice 3' splice sites for each shRNA-seq.

 `star_abs5_motif` stores the motifs of cryptic splice 5' splice sites for each shRNA-seq.

 `star_abs3_motif_canonical` stores the motifs of cryptic splice 3' splice sites near canonical splice sites for each shRNA-seq.

 `star_abs5_motif_canonical` stores the motifs of cryptic splice 5' splice sites near canonical splice sites for each shRNA-seq.

 `star_abs3_motif_exon` stores the motifs of cryptic splice 3' splice sites in exon regions for each shRNA-seq.

 `star_abs5_motif_exon` stores the motifs of cryptic splice 5' splice sites in exon regions for each shRNA-seq.

 `star_abs3_motif_deep_intron` stores the motifs of cryptic splice 3' splice sites in deep introns for each shRNA-seq.

 `star_abs5_motif_deep_intron` stores the motifs of cryptic splice 5' splice sites in deep introns for each shRNA-seq.


### Metaprofiles to show cryptic splice sites' postion (.pdf)
The `pdf` file in each cell line's directory
The meta-profiles of cryptic sites in K562, “3_” indicates cryptic 3' splice sites, while “5_” indicates cryptic 5' splice sites. Three panels in each figure: the first panel is the metaprofile of 5' splice sites, where 0 is the 5' splice site position; the second panel is the metaprofile of 3' splice sites, where 0 is the 3' splice site position; the third panel is around the brand point, where 0 is the branch point position.


### Suggestions and comments are very welcome: limeng49631@aliyun.com



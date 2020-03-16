library(stringr);
library(plyr)
library(dplyr)
library(idr)
setwd("/Users/mengli/Documents/projects/abs/");

data<-read.table("mouse/g_53_sj_all_read_count_mouse.tsv",
                 header = TRUE,as.is = TRUE,sep = "\t",comment.char="$",quote = "$");

data[,"corrected.cryptic.5.ss"]<-data[,"cryptic.5.ss"]/data[,"read_count"];


data_cgr8<-data[data$is_N2A==FALSE,c("rbp","cryptic.5.ss","cryptic.3.ss","corrected.cryptic.3.ss","corrected.cryptic.5.ss")];

data_n2a<-data[data$is_N2A==TRUE,c("rbp","cryptic.5.ss","cryptic.3.ss","corrected.cryptic.3.ss","corrected.cryptic.5.ss")];


colnames(data_cgr8)<-c("rbp","cryptic.5.ss.cgr8","cryptic.3.ss.cgr8","corrected.cryptic.3.ss.cgr8","corrected.cryptic.5.ss.cgr8");

colnames(data_n2a)<-c("rbp","cryptic.5.ss.n2a","cryptic.3.ss.n2a","corrected.cryptic.3.ss.n2a","corrected.cryptic.5.ss.n2a");



data_cgr8[,"rbp"]<-str_sub(data_cgr8[,"rbp"],1,nchar(data_cgr8[,"rbp"])-5);

data_n2a[,"rbp"]<-str_sub(data_n2a[,"rbp"],1,nchar(data_n2a[,"rbp"])-4);


data_cgr8_n2a<-inner_join(data_cgr8,data_n2a,by=c("rbp"="rbp"));



mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.02
idr.out <- est.IDR(data_cgr8_n2a[,c("corrected.cryptic.3.ss.cgr8","corrected.cryptic.3.ss.n2a")], 
                   mu, sigma, rho, p, eps=0.001, max.ite=100);

idr_3ss<-idr.out$IDR

data_cgr8_n2a<-cbind(data_cgr8_n2a,idr_3ss);


mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.002
idr.out <- est.IDR(data_cgr8_n2a[,c("corrected.cryptic.5.ss.cgr8","corrected.cryptic.5.ss.n2a")], 
                   mu, sigma, rho, p, eps=0.001, max.ite=200);
idr_5ss<-idr.out$IDR

data_cgr8_n2a<-cbind(data_cgr8_n2a,idr_5ss);

write.table(data_cgr8_n2a,"mouse/cgr8_n2a_idr.tsv",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


# idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.01, max.ite=20)


# 8245.05
# 7714.325



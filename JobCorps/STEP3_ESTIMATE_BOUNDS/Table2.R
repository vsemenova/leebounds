rm(list=ls())
library(xtable)


###  save your path
name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}
if (name == "Mengsi") {
  my_path<-"C:/Users/mengs/Dropbox/leebounds/"
  # my_path<-"/Users/mac/Dropbox/leebounds/"
}


source(paste0(my_path,"/R/auxiliary.R"))

at_share34<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lee_ATSHARE.csv"))
at_share56<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Full_ATSHARE.csv"))
at_share78<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lasso_ATSHARE.csv"))


estimates2<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lee_estimates.csv"),row.names=1)
estimates3<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Full_estimates.csv"),row.names=1)
estimates4<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lasso_estimates.csv"),row.names=1)

sd2<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lee_CR.csv"),row.names=1)
sd3<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Full_CR.csv"),row.names=1)
sd4<-read.csv(paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/csv/Table1_Lasso_CR.csv"),row.names=1)

estimates<-t(as.matrix(cbind(estimates2, estimates3, estimates4)))
sd<-t(as.matrix(cbind(sd2,sd3,sd4)))

estimates_basic<-t(as.matrix(cbind(estimates2[,1], estimates3[,1], estimates4[,1])))
sd_basic<-t(as.matrix(cbind(sd2[,1],sd3[,1],sd4[,1])))

table4latex<-print_table(estimates,sd,digs=3)
write.table(print(xtable(table4latex, type="latex"),include.rownames =FALSE ),paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/txt/Table2.txt"),append=FALSE)

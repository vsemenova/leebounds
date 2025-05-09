## define a set of covariates
rm(list=ls())

library(sas7bdat)
library(dplyr)
library(tidyverse)
library(feather)

my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/JobCorps/"

baseline<-read.sas7bdat(paste0(my_path,"/Raw_Data/baseline.sas7bdat"))
Lee_data<-read.csv(paste0(my_path,"/Derived_Data/dataLee2009.csv"))
Lee_data[is.na(Lee_data)]<-0 ## ignore NAs

numeric_names<-colnames(baseline)[sapply(baseline,is.numeric)]
factor_names<-colnames(baseline)[sapply(baseline,is.factor)]
true_numeric_names<-c()
# max number of levels for discretization
max_factor<-15
numeric_names<-intersect(numeric_names,colnames(Lee_data))

Lee_data[,numeric_names]<-sapply(Lee_data[,numeric_names],round,0)

for (name in numeric_names) {
  
  x<-length(unique(as.numeric(unlist(as.data.frame(Lee_data[,name]))) ))
  
  if (x<=max_factor & x>=2) {
    factor_names<-c(name,factor_names)
  } else {
    true_numeric_names<-c(name,true_numeric_names)
  }
  
}

factor_names<-intersect(factor_names,colnames(Lee_data))
factor_names<-setdiff(factor_names,"MPRID")

Lee_data[,factor_names]<-sapply(Lee_data[,factor_names],as.factor)
Lee_data[,factor_names][is.na(Lee_data[,factor_names])]<-0
Lee_data_factors_expanded<-model.matrix(~.,Lee_data[,factor_names])
Lee_data_numeric<-Lee_data[,true_numeric_names]
Lee_data_covariates<-cbind(MPRID=Lee_data$MPRID,Lee_data_factors_expanded,Lee_data_numeric)
names(Lee_data_covariates)<-make.unique(names(Lee_data_covariates))
colnames(Lee_data_covariates)[1]<-"MPRID"
selected_names<-setdiff(colnames(Lee_data_covariates),c("(Intercept)","X.Intercept."))
Lee_data_covariates<-Lee_data_covariates[,selected_names]
# too large for GitHub

work_names<-setdiff( (colnames(Lee_data))[startsWith(names(Lee_data),"WORK")], 
                     (colnames(Lee_data))[startsWith(names(Lee_data),"WORKH")])
welf_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"WELF")]
trng_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"TRNG")]
schl_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"SCHL")]
sercr_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"SERCR")]

Lee_data_factors_expanded<-model.matrix(~.,Lee_data[,c(work_names,welf_names,trng_names, schl_names)])
Lee_data_covariates2<-cbind(MPRID=Lee_data$MPRID,Lee_data_covariates,Lee_data_factors_expanded)
names(Lee_data_covariates2)<-make.unique(names(Lee_data_covariates2))
colnames(Lee_data_covariates2)[1]<-"MPRID"
selected_names<-setdiff(colnames(Lee_data_covariates2),c("(Intercept)","X.Intercept."))

### test: the data should contain WORK252, WORK3023 (all those covariates you did not find last time)
write_feather(Lee_data_covariates2,paste0(my_path,"/Derived_Data/dataLee2009covariates.feather"))
## 
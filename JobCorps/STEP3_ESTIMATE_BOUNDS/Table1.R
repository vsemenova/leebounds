## Table 1: Colummns (1) and (2) bounds under unconditional monotonicity

rm(list=ls())
library(tidyverse)
library(dplyr)
library(expm)
library(questionr)

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}
if (name == "Mengsi") {
  my_path<-"C:/Users/mengs/Dropbox/leebounds/"
  # my_path<-"/Users/mac/Dropbox/leebounds/"
}

## load data
Lee_data<-read.csv(paste0(my_path,"JobCorps/Derived_Data/dataLee2009.csv"))

## attach helper functions
source(paste0(my_path,"R/auxiliary.R"))
source(paste0(my_path,"R/leebounds.R"))

#sink(paste0(my_path,"JobCorps/Logs/STEP3_Bounds/Tab1_Cols12.log"))

### number of bootstrap repetitions
print("Number of Bootstrap Repetitions")
### other parameters
N_rep=1000

print ("Number of Boot rep ...")
print (N_rep)
ci_alpha=0.05

print ("Conf Level ...")
print (ci_alpha)

## main arguments 
week<-90
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")

### create dataset for week = week
leedata_cov<-prepare_leedata(week=week,Lee_data,covariate_data = Lee_data[, baseline_varnames   ])

##  basic Lee bound (Column 1)
res_basic<-basic_lee_bound(leedata_cov,treat_helps=TRUE)
bounds_basic<-GetBounds(res_basic)
print("Basic Lee Bounds")
print(bounds_basic)
print("Always-takers' share")
print(res_basic$p0)

bounds_bb<-main_bb(function_name=basic_lee_bound,
                   mydata=leedata_cov[,c("treat","selection","outcome","weights")],
                   N_rep=N_rep,treat_helps=TRUE)

## confidence region for identified set
print ("Confidence region for Basic Bounds")

basic_CR<-compute_confidence_region(bounds_bb,bounds_basic, ci_alpha=ci_alpha)

print(basic_CR)

### discrete Lee bound (Column 2)
### construct predicted wage potential covariates
### STEP 2.1: create wage potential covariate
lm_model<-lm(outcome~.,data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,c("outcome",baseline_varnames)])
predicted_wage <-  cbind(rep(1,dim(leedata_cov)[1]), as.matrix(leedata_cov[,baseline_varnames]))%*%lm_model$coefficients

## cutoffs for predicted wage as in Lee (2009)  
intervals <- c(-Inf, log(6.75), log(7), log(7.5), log(8.5), Inf)

### STEP 2.2. discretize wage potential into 5 groups
wage_groups <- cut(predicted_wage, breaks = intervals, labels = c("Low", "Low-Mid", "Middle", "Mid-High", "High"), include.lowest = TRUE)
### take group as the indicator of predicted wage group
leedata_cov$group<- wage_groups


### STEP 2.3. Estimate bounds, enforcing S(1) \geq S(0) for each group, i.e. treat_helps=TRUE
res_discrete<-discrete_bound(leedata_cov,treat_helps=TRUE)
bounds_discrete<-GetBounds(res_discrete)

print("Discrete Lee Bounds")
print(bounds_discrete)
bounds_bb<-main_bb(function_name=discrete_bound,
                   mydata=leedata_cov[,c("treat","selection","outcome","weights","group")],N_rep=N_rep,treat_helps=TRUE)
## confidence region for identified set
discrete_CR<-compute_confidence_region(bounds_bb,bounds_discrete, ci_alpha=ci_alpha)
print("Confidence Region for Discrete Lee Bounds")
print(discrete_CR)

### save results
results<-matrix(0,2,2)
results[,1]<-round(bounds_basic,3)
results[,2]<-round(bounds_discrete,3)
colnames(results)<-c("basic", "discrete")
rownames(results)<-c("lower","upper")

results_CR<-matrix(0,2,2)
results_CR[,1]<-round(basic_CR,3)
results_CR[,2]<-round(discrete_CR,3)
colnames(results_CR)<-c("basic", "discrete")
rownames(results_CR)<-c("lower","upper")

closeAllConnections()


estimates<-t(as.matrix(cbind(results)))
sd<-t(as.matrix(cbind(results_CR)))
table4latex<-print_table(estimates,sd,digs=3)
write.table(print(xtable(table4latex, type="latex"),include.rownames =FALSE ),
            paste0(my_path,"JobCorps/Tables/STEP3_Bounds/",name,"/txt/Table1.txt"),append=FALSE)


# Print LaTeX table



#write.csv(results, paste0(my_path,"JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Col12_estimates.csv"))
#write.csv(results_CR, paste0(my_path,"JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Col12_CR.csv"))
#write.csv(res_basic$p0, paste0(my_path,"JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Col1_ATSHARE.csv"))
#write.csv(res_discrete$p0, paste0(my_path,"JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Col2_ATSHARE.csv"))
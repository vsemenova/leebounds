rm(list=ls())
library(tidyverse)
library(dplyr)
library(hdm)
library(expm)

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
} else {
  my_path<-"CHANGE_ME"
}

Lee_data<-read.csv(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009.csv"))

source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/first_stage_functions.R"))
source(paste0(my_path,"/R/second_stage_functions.R"))

## main arguments 
week<-90
thresh<-0.1
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
leedata_cov<-prepare_leedata(week=week,Lee_data,covariate_data = Lee_data[,baseline_varnames])


### other parameters
N_rep=100

print("Number of Boot rep ...")
print(N_rep)
ci_alpha=0.05

print("Conf Level ...")
print(ci_alpha)

quantile_grid_size=0.01
min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
sample_size<-dim(leedata_cov)[1]
rhoN = 1*(sample_size)^{-1/4} * (log ( (sample_size)))^{-1}
p0 =1-rhoN
p1=1+rhoN

### Column (5)

form_selection<-as.formula(paste0("selection~(treat)*(", paste0(baseline_varnames,collapse="+"),")"))
glm.fit<-estimate_selection(form=form_selection,leedata=leedata_cov,selection_function_name = "glm")


###  basic generalized bound
s.hat<-as.data.frame(predict_selection(glm.fit, leedata=leedata_cov))
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
leedata_cov$p.0.star<-p.0.star
leedata_cov$s.0.hat<-s.hat$s.0.hat
leedata_cov$s.1.hat<-s.hat$s.1.hat

inds_not_boundary<-!((p.0.star>1-thresh) & (p.0.star<1+thresh))
bounds_basic<-GetBounds(basic_generalized_bound(leedata=leedata_cov[inds_not_boundary,]))
print ("Basic Generalized Bounds ...")
print(bounds_basic)

bounds_bb<-main_bb(function_name=basic_generalized_bound,
                   mydata=leedata_cov[inds_not_boundary,c("treat","selection","outcome","weights","p.0.star","s.0.hat","s.1.hat","prop1","prop0")],
                   N_rep=N_rep)
## confidence region for identified set
print ("Confidence Region for Basic Generalized Bounds ...")
basic_CR<-compute_confidence_region(bounds_bb,bounds_basic, ci_alpha=0.05 )
print(basic_CR)

res_tight<-ortho_leebounds(leedata=leedata_cov,
                           s.hat=s.hat,
                           quantile_grid_size = quantile_grid_size,
                           variables_for_outcome=baseline_varnames,
                           min_wage=min_wage,
                           max_wage=max_wage,
                           ortho=TRUE,
                           p0=1-rhoN, p1=1+rhoN,outcome_function_name="rq")

bounds_tight<-GetBounds(res_tight)

print ("Estimated bounds")
print(bounds_tight)
print ("Estimated always-takers' share")
print(res_tight$denom)

bounds_bb<-main_bb(res_tight$leedata,N_rep=N_rep,function_name=second_stage_wrapper,ortho=TRUE,p0=1-rhoN, p1=1+rhoN)
tight_CR<-compute_confidence_region(ATE_boot=bounds_bb,ATE_est= bounds_tight,ci_alpha=ci_alpha)

### save results
results<-matrix(0,2,2)
results[,1]<-round(bounds_basic,3)
results[,2]<-round(bounds_tight,3)
colnames(results)<-c("basic", "tight")
rownames(results)<-c("lower","upper")

print("Estimated Bounds")
print(results)

results_CR<-matrix(0,2,2)
results_CR[,1]<-round(basic_CR,3)
results_CR[,2]<-round(tight_CR,3)
colnames(results_CR)<-c("basic", "tight")
rownames(results_CR)<-c("lower","upper")

print("Estimated CR")
print(results_CR)

write.csv(results, paste0(my_path,"/JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Lee_estimates.csv"))
write.csv(results_CR, paste0(my_path,"/JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Lee_CR.csv"))
write.csv(print(res_tight$denom), paste0(my_path,"/JobCorps/Tables/STEP3_Bounds/", name, "/csv/Table1_Lee_ATSHARE.csv"))


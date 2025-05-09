rm(list=ls())

library(feather)

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}
if (name == "Mengsi") {
  # my_path<-"C:/Users/mengs/Dropbox/leebounds/"
  my_path<-"/Users/mac/Dropbox/leebounds"
}

Lee_data<-read.csv(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)

source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/first_stage_functions.R"))
source(paste0(my_path,"/R/second_stage_functions.R"))

TREAT<-Lee_data$TREATMNT.y
# covariates selected by David Lee (2009)
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]
#Lee_data_all_covariates<-cbind(Lee_data[,baseline_varnames],Lee_data_all_covariates)

weeks<-1:208
emp_rate_Xhelp_treat<-rep(0,length(weeks))
emp_rate_Xhelp_control<-rep(0,length(weeks))
emp_rate_Xhurt_treat<-rep(0,length(weeks))
emp_rate_Xhurt_control<-rep(0,length(weeks))

emp_rate_X_treat<-rep(0,length(weeks))
emp_rate_X_control<-rep(0,length(weeks))

estimated.p.0.hat.nonmonotone<-matrix(0,dim(Lee_data_covariates),length(weeks))
emp_rate_diff<-rep(0,length(weeks))

p.0.star<-matrix(0,dim(Lee_data_covariates),length(weeks))

form_nonmonotone<-as.formula(paste0("selection~(treat)*(", paste0(baseline_varnames,collapse="+"),")"))

weeks<-1:208
s.hat.glm<-array(0,c(dim(Lee_data_covariates)[1],2,length(weeks)))
fraction_treat_helps<-rep(0,length(weeks))
s_min<-rep(0,length(weeks))
s_max<-rep(0,length(weeks))

for (i in 10:length(weeks)) {
  week<-weeks[i]
  print(paste0("Week ",i))
  
  leedata_cov<-prepare_leedata(week,Lee_data,covariate_data=Lee_data_covariates)
  
  glm.fit<-estimate_selection(form_nonmonotone,leedata=leedata_cov, selection_function_name = "glm")
  
  res<-predict_selection(glm.fit,leedata_cov)
  s.hat = data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  
  p.0.star[,i]<-s.hat$s.0.hat/s.hat$s.1.hat
  
  tau.x<-s.hat$s.1.hat - s.hat$s.0.hat
  
  inds_helps<-p.0.star[,i]<1
  inds_hurts<-!(inds_helps)
  
  emp_rate_Xhelp_treat[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==1 & inds_helps ],
                                         weights=leedata_cov$weights[leedata_cov$treat==1 & inds_helps ])
  emp_rate_Xhelp_control[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==0 & inds_helps ],
                                           weights=leedata_cov$weights[leedata_cov$treat==0 & inds_helps ])
  
  emp_rate_Xhurt_treat[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==1 & inds_hurts ],
                                         weights=leedata_cov$weights[leedata_cov$treat==1 & inds_hurts ])
  
  emp_rate_Xhurt_control[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==0 & inds_hurts ],
                                           weights=leedata_cov$weights[leedata_cov$treat==0 & inds_hurts ])
  s.hat.glm[,,i]<-as.matrix(s.hat,ncol=2)
  
  
  ## Figure 1: average employment rate
  emp_rate_X_treat[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==1]==1,
                                     weights=leedata_cov$weights[leedata_cov$treat==1])
  
  emp_rate_X_control[i]<-weighted.mean(leedata_cov$selection[leedata_cov$treat==0]==1,
                                       weights=leedata_cov$weights[leedata_cov$treat==0])
  
  emp_rate_diff[i]<-(emp_rate_X_treat[i]- emp_rate_X_control[i])
  
  
  ## Figure 2: 
  fraction_treat_helps[i]<-weighted.mean(p.0.star[,i]<=1,weights=Lee_data$DSGN_WGT.y)
  
  ## Figure 6
  
  leedata_cov$s.0.hat<-s.hat$s.0.hat
  leedata_cov$s.1.hat<-s.hat$s.1.hat
}

Figure1_dataset=data.frame(weeks=c(weeks,weeks),
                           delta=c((emp_rate_Xhelp_treat-emp_rate_Xhelp_control),
                                   (emp_rate_Xhurt_treat-emp_rate_Xhurt_control)),
                           group=c(rep(1,length(weeks)),rep(0,length(weeks))),
                           emp_rate_diff= c(emp_rate_diff,emp_rate_diff))

write.csv(Figure1_dataset,paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/",name,"/Figure1.csv"))

selected_weeks<-weeks
Figure2_dataset<-data.frame(weeks=weeks, fraction=fraction_treat_helps)
write.csv(Figure2_dataset,paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/",name,"/Figure2.csv"))

rm(list=ls())
#install.packages("feather")
library(feather)
library(expm)

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds"
}
if (name == "Mengsi") {
  my_path<-"C:/Users/mengs/Dropbox/leebounds"
}

source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/first_stage_functions.R"))
source(paste0(my_path,"/R/second_stage_functions.R"))

print ("Loading data ...")
Lee_data<-read.csv(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009.csv"))

Lee_data<-as.data.frame(Lee_data)
Lee_data[is.na(Lee_data)]<-0
### FIGURE 1 and 2###
weeks<-1:208
quantile_grid_size=0.01

baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR"  )

### FIGURE3
bounds_non_treated_wage<-matrix(0,2,length(weeks))




for (i in 1:length(weeks)) {
  print(i)
  week<-weeks[i]
  leedata_cov<-prepare_leedata(week,Lee_data,covariate_data=Lee_data[,baseline_varnames])
 
  
  ## analysis 
  glm.fit<-estimate_selection(leedata= leedata_cov,variables_for_selection = baseline_varnames, selection_function_name = "glm")
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  
 
  min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
  max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  
  leebounds_ortho_result<-ortho_leebounds(leedata=leedata_cov,
                                          s.hat=s.hat,
                                          function_ss=numerator_control_wage_2ndstage,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=baseline_varnames,
                                          min_wage=min_wage,
                                          max_wage=max_wage,
                                          ortho=TRUE,
                                          p0=1, p1=1,outcome_function_name="rq", function_ss_name ="control wage")
  
  bounds_non_treated_wage[,i]<-GetBounds(leebounds_ortho_result)
 
}


Figure3_dataset<-data.frame(weeks=c(weeks,
                                    weeks),
                            bound=c(bounds_non_treated_wage[1,],bounds_non_treated_wage[2,]),
                            group=c(rep(0,length(weeks)),rep(1,length(weeks))))
write.csv(Figure3_dataset,paste0(my_path,"/JobCorps/Tables/STEP4_Misc/", name, "/Figure3.csv"))
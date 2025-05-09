rm(list=ls())

library(feather)
library(tidyverse)
library(radiant.data)

### for parallel computing, replacable by standard for loop
library(doParallel)
library(dplyr)
library(foreach)

cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}


source(paste0(my_path,"/JobCorps/STEP2_TEST_MONOTONICITY/utils_for_test.R"))

Lee_data<-read.csv(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009covariates.feather"))

baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")

Lee_data_all_covariates<-cbind(Lee_data_all_covariates,Lee_data[,baseline_varnames])
Lee_data_all_covariates[is.na(Lee_data_all_covariates)]<-0
Lee_data_all_covariates$EARN_YR_quant<-Lee_data_all_covariates$EARN_YR>720 & Lee_data_all_covariates$EARN_YR<3315
Lee_data_all_covariates$EARN_YR_perc<-Lee_data_all_covariates$EARN_YR>3315 & Lee_data_all_covariates$EARN_YR<7298

alpha<-0.05
crit.val.05<-qnorm(1-alpha/2)/sqrt(1-(qnorm(1-alpha/2))^2/dim(Lee_data)[1])

alpha<-0.01
crit.val.01<-qnorm(1-alpha/2)/sqrt(1-(qnorm(1-alpha/2))^2/dim(Lee_data)[1])


# week=60:89
cov_names<-c("PERS_INC3","EARN_YR_quant","MOS_AFDC8")
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)
mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,
                        Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=30]<-2
mygroup<-ungroup(mygroup)
mygroup$group[mygroup$group %in% c(0,1,2)]<-0

group_des_M0 <- unique(mygroup[, -c(2, 3)])

myres_M0=foreach(week=60:89, .combine = 'cbind',.packages = c("tidyverse","radiant.data"))  %dopar% {
  ## constructed data for week
  res<-test_wrapper(week,mygroup)
  res
}

# week=90:116
## constructed data for week
cov_names<-c("TYPEWORR5","EARN_YR_perc","R_HOME1","MARRCAT11","R_COMM1","WELF_KID4")
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,
                        Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=30]<-1
mygroup<-ungroup(mygroup)
mygroup$group[mygroup$group !=7]<-0

group_des_M1 <- unique(mygroup[, -c(2, 3)])

myres_M1=foreach(week=90:116, .combine = 'cbind',.packages = c("tidyverse","radiant.data"))  %dopar% {
  res<-test_wrapper(week,mygroup)
}

# week=117:152
## constructed data for week
cov_names<-c("EARN_YR_quant","R_HOME1","TYPEWORR5","R_COMM1","FRQ_POT3","DRG_SUMP2","IMP_PRO1","MARRCAT11")

group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,
                        Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=30]<-1
mygroup$group[mygroup$group !=16]<-0

group_des_M2 <- unique(mygroup[, -c(2, 3)])

myres_M2=foreach(week=117:152, .combine = 'cbind',.packages = c("tidyverse","radiant.data"))  %dopar% {
  res<-test_wrapper(week,mygroup)
  res
}

# week=153:186
## constructed data for week
cov_names<-c("EARN_YR_quant","R_HOME1","TYPEWORR5","R_COMM1","FRQ_POT3","DRG_SUMP2","REASED_R4","IMP_PRO1","MARRCAT11")

group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,
                        Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=30]<-1
mygroup$group[mygroup$group !=17]<-0 

group_des_M3 <- unique(mygroup[, -c(2, 3)])

myres_M3=foreach(week=153:186, .combine = 'cbind',.packages = c("tidyverse","radiant.data"))  %dopar% {
  res<-test_wrapper(week,mygroup)
  res
}

# week=187:208
## constructed data for week
cov_names<-c("TYPEWORR5","EARN_YR_perc","R_HOME1","MARRCAT11","R_COMM1","WELF_KID4")
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,
                        Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=30]<-1
mygroup<-ungroup(mygroup)
mygroup$group[mygroup$group !=7]<-0

group_des_M4 <- unique(mygroup[, -c(2, 3)])

myres_M4=foreach(week=187:208, .combine = 'cbind',.packages = c("tidyverse","radiant.data"))  %dopar% {
  res<-test_wrapper(week,mygroup)
  res
}


test_result<-data.frame(weeks=c(1:208),tstat=rep(1,208),pvalue=rep(1,208))
test_result$tstat[60:89]<-as.numeric(apply(myres_M0[1:2, ], 2, max))
test_result$tstat[90:116]<-as.numeric(apply(myres_M1[1:2, ], 2, max))
test_result$tstat[117:152]<-as.numeric(apply(myres_M2[1:2, ], 2, max))
test_result$tstat[153:186]<-as.numeric(apply(myres_M3[1:2, ], 2, max))
test_result$tstat[187:208]<-as.numeric(apply(myres_M4[1:2, ], 2, max))


find_my_pvalue<-function(x) {
  if (x>=crit.val.01) {
    p<-0.01
  } else if (x>=crit.val.05){
    p<-0.05
  } else {
    p<-1
  }
  return(p)
}
test_result$pvalue<-sapply(test_result$tstat,find_my_pvalue)
write.csv(test_result,paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/",name,"/test_result.csv"))


### summarize Table D.11
sink(paste0(my_path,"/JobCorps/Logs/step2_t_stats.log"))
# Print the means
cat("week 60-89:", mean(test_result$tstat[60:89]), "\n")
cat("week 90-116:", mean(test_result$tstat[90:116]), "\n")
cat("week 117-152:", mean(test_result$tstat[117:152]), "\n")
cat("week 153-186:", mean(test_result$tstat[153:186]), "\n")
cat("week 187-208:", mean(test_result$tstat[187:208]), "\n")


closeAllConnections()



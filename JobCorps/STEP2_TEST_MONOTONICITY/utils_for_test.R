compute_tstat_by_group<-function(leedata,...) {
  
  # normalization by 1/2 does not matter for t-ratio
  # by weighting, we equalize prob of being treated = prob of being untreated = 0.5
  X_inds<- ((leedata$selection ==1 ) * (leedata$treat==1) - 
              (leedata$selection ==1 ) * (leedata$treat==0))/0.5
  leedata<-as.data.frame(leedata)
  
  tstat<-sqrt(dim(leedata)[1])* weighted.mean(  X_inds,w=leedata$weights)/weighted.sd(  X_inds,wt=leedata$weights)
  
  return(tstat)
}


test_wrapper<-function(week,mygroup){
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  if (week>=90) {
    leedata_cov<-data.frame(treat=1-Lee_data$TREATMNT.y,
                            selection = logwage_week>0,
                            MPRID=Lee_data$MPRID,
                            weights=Lee_data$DSGN_WGT.y) %>%
      inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
    
  } else {
    leedata_cov<-data.frame(treat=Lee_data$TREATMNT.y,
                            selection = logwage_week>0,
                            MPRID=Lee_data$MPRID,
                            weights=Lee_data$DSGN_WGT.y) %>%
      inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
    
  }
  ## created grouped data
  grouped_data<-group_by(leedata_cov,group)
  
  ### compute test statistic
  grouped_data<-group_by(grouped_data,group)
  grouped_data_tstat<-grouped_data
  tstat<-unlist(group_map(grouped_data_tstat,compute_tstat_by_group,keep=TRUE))
  max.t.stat<-max(tstat)
  index_of_max_tstat <- which.max(tstat)
  
  result <- c(tstat, index_of_max_tstat)
  
  return(result)
}



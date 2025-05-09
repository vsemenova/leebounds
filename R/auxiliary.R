### create a cross-sectional data (D,S,X,SY) from Lee_data and covariate_data
### attach weights
prepare_leedata<-function(week,Lee_data,covariate_data,prop1=NULL,prop0=NULL,...){
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome=logwage_week))
  
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome = as.numeric(Lee_data[,earn_name]>0)*logwage_week)
  
  if (!is.null(covariate_data)) {
    leedata_cov<-cbind(leedata_cov, covariate_data)
  }
  
  leedata_cov[is.na(leedata_cov)]<-0
  
  leedata_cov$weights<-Lee_data$DSGN_WGT.y
  
  if (is.null(prop1)) {
    leedata_cov$prop1<-1/2
    leedata_cov$prop0<-1/2
  } else {
    leedata_cov$prop1<-prop1
    leedata_cov$prop0<-prop0
  }
  
  return(leedata_cov)
}

## method to get bounds
main_bb<-function(mydata,function_name,N_rep=10,...) {
  
  ATE_bb<-matrix(0,N_rep,2)
  sample_size<-dim(mydata)[1]
  for (b in 1:N_rep) {
    set.seed(b)
    #print(b)
    inds<-sample(1:sample_size,sample_size,replace=TRUE)
    mydatab<-mydata[inds,]
    resultb = try(function_name(mydatab,...))
    # resultb = try(function_name(mydatab))
    ATE_bb[b,]<-GetBounds(resultb)
  }
  
  return(ATE_bb)
}


compute_confidence_region<-function(ATE_boot,ATE_est,ci_alpha=0.05,tol=1e-5) {
  Omega.hat<-matrix(0,2,2)
  if (sum(is.na(ATE_boot))+sum(is.na(ATE_est))>0) {
    return(c(lower_bound = NA, upper_bound=NA))
  }
  ATE_boot_centered<-matrix(0,dim(ATE_boot)[1],2)
  ## Centered draws of lower bound
  ATE_boot_centered[,1]<-ATE_boot[,1]-ATE_est[1]
  ## Centered draws of upper bound
  ATE_boot_centered[,2]<-ATE_boot[,2]-ATE_est[2]
  
  Omega.hat[1,1]<-var( ATE_boot_centered[,1])
  Omega.hat[2,2]<-var(ATE_boot_centered[,2])
  Omega.hat[1,2]<-cov(ATE_boot_centered[,1],ATE_boot_centered[,2])
  Omega.hat[2,1]<-Omega.hat[1,2]
  
  crit.val<-sqrtm(Omega.hat)%*% c(-qnorm(sqrt(1-ci_alpha)), qnorm(sqrt(1-ci_alpha)) ) 
  if (max(abs(Im(sqrtm(Omega.hat))))>tol) {
    stop ("Non-trivial imaginary part!")
  } else {
    crit.val<-sapply( crit.val,Re)
    
  }
  lower_bound<-ATE_est[1]+ crit.val[1]
  upper_bound<-ATE_est[2] +crit.val[2]
  return(c(lower_bound = lower_bound, upper_bound=upper_bound))
}


GetBounds<-function(x) {
  return(c(x$lower_bound,x$upper_bound))
}

GetFraction<-function(x) {
  return (x$fraction)
}

GetTrimmedBounds<-function(x) {
  return(c(x$trimmed_mean_lower,x$trimmed_mean_upper))
}

GetThresh<-function(x) {
  return(x$p0)
}

GetOdds<-function(x) {
  return(x$odds)
}

Getyp0<-function(x) {
  return(x$yp0)
}

Gety1p0<-function(x) {
  return(x$y1p0)
}

Gets0<-function(x) {
  return(x$s0)
}

Gets1<-function(x) {
  return(x$s1)
}

Getprop0<-function(x) {
  return(x$prop0)
}

Getprop1<-function(x) {
  return(x$prop1)
}
standardize<-function(vec) {
  if (sd(vec)>0) {
    vec<-(vec-mean(vec))/sd(vec)
  } 
  return(vec)
}



print_table<-function(estimates,sd,im=NULL,digs=3) {
  
  estimates<-apply(estimates,2,round,digs)
  sd<-apply(sd,2,round,digs)
  
  M<-matrix(NA,2,dim(estimates)[1])
  

    for (k in 1:(dim(estimates)[1])) {
      print(k)
      M[1,k]<-paste0("[", estimates[k,1],", " ,estimates[k,2],"]")
      M[2,k]<-paste0("(", sd[k,1],", " ,sd[k,2],")")
      
    }
  
  
  
  return(M)
}

### first stage

orthogonal_correction<-function(leedata,treat_helps,...) {
  d<-leedata$treat
  s<-leedata$selection
  sy<-s*leedata$outcome
  
  ## compute second stage estimate based on the first stage
  ## args: first-stage estimate
  y.p0.hat<-leedata$y.p0.hat
  y.1.p0.hat<-leedata$y.1.p0.hat
  s.0.hat<-leedata$s.0.hat
  s.1.hat<-leedata$s.1.hat
  
  prop1<-leedata$prop1
  prop0<-leedata$prop0
  
  if (treat_helps) {
    
    p.0.hat<-s.0.hat/s.1.hat
    p.0.hat<-sapply(p.0.hat,min,0.99999)
    
    gamma1x<-y.1.p0.hat
    gamma2x<- (-1)*(y.1.p0.hat)*p.0.hat
    gamma3x<-(y.1.p0.hat)
    
    alpha1x<- (1-d)/prop0*(s-s.0.hat)
    alpha2x<- d/prop1*(s-s.1.hat)
    alpha3x<- d*s/prop1*(as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat))
    
    gamma4x<-y.p0.hat
    gamma5x<-(-1)*(y.p0.hat)*p.0.hat
    gamma6x<-(-1)*y.p0.hat
    
    alpha4x<- alpha1x
    alpha5x<- alpha2x
    alpha6x<- d*s/prop1*(as.numeric(sy<=y.p0.hat) - p.0.hat)
    
    A1<-gamma1x*alpha1x
    A2<-gamma2x*alpha2x
    A3<-gamma3x*alpha3x
    
    A4<-gamma4x*alpha4x
    A5<-gamma5x*alpha5x
    A6<-gamma6x*alpha6x
    
    ## correction term in Definition 4.1 , equation (4.6)
    upper_trim_correction<-(A1+A2+A3)
    ## correction term in Definition 4.3,  equation (4.11)
    lower_trim_correction<-(A4+A5+A6)
    
    
  } else {
    
    # this expression is 1/p(X) in equation (4.9)
    p.0.hat<-s.1.hat/s.0.hat
    p.0.hat<-sapply(p.0.hat,min,0.99999)
    
    ## Riesz representers functions
    gamma1x<-(y.p0.hat)*p.0.hat
    gamma2x<-(-1)*y.p0.hat
    gamma3x<-(y.p0.hat)
    
    # residuals
    alpha1x<-(1-d)/prop0*(s-s.0.hat)
    alpha2x<-d/prop1*(s-s.1.hat)
    alpha3x<-(1-d)/prop0*s*(as.numeric(sy<=y.p0.hat)- p.0.hat)
    
    ## correction terms
    A1<-gamma1x*alpha1x
    A2<-gamma2x*alpha2x
    A3<-gamma3x*alpha3x
    
    ## correction term in Definition 4.2 , equation (4.9)
    lower_trim_correction<- (A1+A2+A3)
    
    
    
    
    gamma4x<-y.1.p0.hat*p.0.hat
    gamma5x<-(-1)*y.1.p0.hat
    gamma6x<-(-1)*y.1.p0.hat
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    alpha6x<-(1-d)/prop0*s*(as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat))
    
    
    
    A4<-gamma4x*alpha4x
    A5<-gamma5x*alpha5x
    A6<-gamma6x*alpha6x
    
    ## correction term in Definition 4.3,  equation (4.12)
    upper_trim_correction<-(A4+A5+A6)
    
    
  }
  
  
  return(list(lower_trim_correction=lower_trim_correction,
              upper_trim_correction=upper_trim_correction))
}


round_p0<-function(p.0.star,thresh=0.94,...) {
  

  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  if (mean(inds_helps)>=thresh) {
    p.0.star<-sapply(p.0.star,min,0.9999)
  } 
  
  if (mean(inds_hurts)>=thresh) {
    print("cutoff")
    p.0.star<-sapply(p.0.star,max,1.0001)
  } 
  return(p.0.star)
}


### wrapper functions to estimate first and second stage
trimmed_moment<-function(leedata,treat_helps,...) {
  ## args: data
  d<-leedata$treat
  s<-leedata$selection
  sy<-s*leedata$outcome
  prop1<-leedata$prop1
  prop0<-leedata$prop0

  ## args: first-stage estimate
  y.p0.hat<-leedata$y.p0.hat
  y.1.p0.hat<-leedata$y.1.p0.hat
  
  if (treat_helps) {
    trimmed_mean_upper<-(d*s*sy*(sy>=y.1.p0.hat))/prop1
    trimmed_mean_lower<-(d*s*sy*(sy<=y.p0.hat))/prop1
    
  } else {
    trimmed_mean_lower<-((1-d)*s*sy*(sy<=y.p0.hat))/prop0
    trimmed_mean_upper<-((1-d)*s*sy*(sy>=y.1.p0.hat))/prop0
  }
  return(list(trimmed_mean_lower=trimmed_mean_lower,
              trimmed_mean_upper=trimmed_mean_upper))
}


numerator_2ndstage<-function(leedata,treat_helps, ortho=TRUE,...) {
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-s*leedata$outcome
  prop1<-leedata$prop1
  prop0<-leedata$prop0
  
  if (ortho==TRUE) {
    correction<-orthogonal_correction(leedata=leedata,treat_helps=treat_helps,...)
  } else {
    correction<-list(lower_trim_correction=0,upper_trim_correction=0)
  }
  
  trimmed_moment_res<-trimmed_moment(treat_helps=treat_helps,leedata=leedata,...)
  trimmed_mean_lower<-trimmed_moment_res$trimmed_mean_lower
  trimmed_mean_upper<-trimmed_moment_res$trimmed_mean_upper
  
  if (treat_helps) {
    moment_upper<-trimmed_mean_upper + correction$upper_trim_correction - sy*(1-d)/prop0 
    moment_lower<-trimmed_mean_lower + correction$lower_trim_correction - sy*(1-d)/prop0
  }
  else {
    moment_upper<-sy*(d)/prop1 - trimmed_mean_lower + correction$lower_trim_correction
    moment_lower<-sy*(d)/prop1 - trimmed_mean_upper + correction$upper_trim_correction
  }
  
  return (list(lower_bound=moment_lower,
               upper_bound=moment_upper))
}

always_takers_share<-function(leedata, ortho_d=TRUE,...) {
  
  sample_size<-dim(leedata)[1]

  d<-leedata$treat
  s<-leedata$selection
  
  
  if (!is.null(leedata$s.hat)) {
    
    s.hat<-leedata$s.hat
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
  } else {
    s.0.hat<-leedata$s.0.hat
    s.1.hat<-leedata$s.1.hat
  }
 
 
  prop1<-leedata$prop1
  prop0<-leedata$prop0
  
  weights<-leedata$weights
  
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  }
  
  p.0.star<-leedata$p.0.star

  inds_helps<-p.0.star<1
  inds_hurts<-!inds_helps
  
  moment_denom<-rep(0,length(d))
  
  if (ortho_d) {
    moment_denom[inds_helps]<-((1-d)*(s-s.0.hat)/prop0 + s.0.hat)[inds_helps]
    moment_denom[inds_hurts]<-(d*(s-s.1.hat)/prop1 + s.1.hat)[inds_hurts]
  } else {
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    moment_denom<-apply(s.hat,1,min)
  }
    
  #denom<-weighted.mean(moment_denom,weights)
  denom<-mean(moment_denom)
  return(denom)
}


second_stage_wrapper<-function(leedata,tol=1e-5,p0=1,p1=1,function_ss_name=NULL,...) {
  
  if (is.null(function_ss_name)) {
    function_ss<-numerator_2ndstage
    function_ss_name<-"effect"
  } else {
    function_ss<-numerator_control_wage_2ndstage
  }
  
  sample_size<-dim(leedata)[1]
  p.0.star<-leedata$p.0.star
  prop1<-leedata$prop1
  prop0<-leedata$prop0
  
  inds_helps<-(p.0.star<=p0)
  inds_hurts<-(p.0.star>p1)
  inds_boundary<- !(inds_helps + inds_hurts)
  
  if (sum(inds_helps)>0){
    res_helps<-function_ss(leedata=leedata[inds_helps,],treat_helps = TRUE, ...)
   }
  else {
    res_helps<-NULL
    estimated_bounds_helps<-c(0,0)
  }
  
  if (sum(inds_hurts)>0){
    res_hurts<-function_ss(leedata=leedata[inds_hurts,], treat_helps = FALSE,...)
  } else {
    res_hurts<-NULL
    estimated_bounds_hurts<-c(0,0)
  }

  if (sum(inds_boundary)>0) {
    res_boundary<-((leedata$treat*leedata$selection*leedata$outcome)/prop1 - ((1-leedata$treat)*leedata$selection*leedata$outcome)/prop0)[inds_boundary]
  } else {
    res_boundary<-0
  }
  
  denom<-always_takers_share(leedata=leedata, ...)
  
  bounds<-c(0,0)
  
  ### weighted mean
  moment_u<-rep(NA,sample_size)
  moment_u[inds_helps]<-res_helps$upper_bound
  moment_u[inds_hurts]<-res_hurts$upper_bound
  
  
  moment_l<-rep(NA,sample_size)
  moment_l[inds_helps]<-res_helps$lower_bound
  moment_l[inds_hurts]<-res_hurts$lower_bound
  
  if (function_ss_name=="effect") {
    moment_u[inds_boundary]<-res_boundary
    moment_l[inds_boundary]<-res_boundary
  }


  if (function_ss_name=="control wage") {
    moment_l[inds_boundary]<-moment_u[inds_boundary]<-(((1-leedata$treat)*leedata$selection*leedata$outcome)/prop0)[inds_boundary]
  }
  
  if (is.null(leedata$weights)) {
    leedata$weights<-1
  }
  
  bounds[1]<-weighted.mean(moment_l,leedata$weights)/denom
  bounds[2]<-weighted.mean(moment_u,leedata$weights)/denom
   
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              denom=denom,
              leedata=leedata
              ))
}


ortho_leebounds<-function(leedata,s.hat=NULL,...) {
  
  ## estimate first stage selection and quantile fitted values
  leedata_fs<-first_stage_wrapper(leedata,s.hat=s.hat,...)
  
  ## calculate bounds
  res<-second_stage_wrapper(leedata=leedata_fs,...)
  return (res)
  
}



### AUXILIARY
numerator_control_wage_2ndstage<-function(leedata,treat_helps, ortho=TRUE,...) {
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-s*leedata$outcome
  prop0<-leedata$prop0
  prop1<-leedata$prop1
  if (ortho==TRUE) {
    correction<-orthogonal_correction(leedata=leedata,treat_helps=treat_helps,...)
  } else {
    correction<-list(lower_trim_correction=0,upper_trim_correction=0)
  }
  
  
  
  if (treat_helps) {
    moment_upper<- sy*(1-d)/prop0 
    moment_lower<- sy*(1-d)/prop0
  }
  else {
    trimmed_moment_res<-trimmed_moment(treat_helps=treat_helps,leedata=leedata,...)
    moment_upper<- trimmed_moment_res$trimmed_mean_upper - correction$upper_trim_correction 
    moment_lower<- trimmed_moment_res$trimmed_mean_lower - correction$lower_trim_correction
  }
  
  return (list(lower_bound=moment_lower,
               upper_bound=moment_upper))
}


frechet_lower_bound<-function(leedata,s.hat,...) {
  
  d<-leedata$treat
  s<-leedata$selection
  
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  
  moment_denom<-sapply(s.hat$s.0.hat+s.hat$s.1.hat-1,max,0)
  
  denom<-mean(moment_denom)
  return(denom)
}

basic_generalized_bound<-function(leedata,ortho_d=TRUE,...) {
  
  d<-leedata$treat
  
  if ("weights" %in% colnames(leedata)) {
    weights<-leedata$weights
  } else {
    weights<-rep(1,length(d))
  }
  p.0.star<-leedata$p.0.star
  
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  if (sum(inds_helps)>0) {
    res_helps <-basic_lee_bound(leedata[inds_helps,],...)
    bounds_helps<-GetBounds(res_helps)
    s_helps<-always_takers_share(leedata[inds_helps,],...)
  } else {
    bounds_helps<-c(0,0)
    s_helps<-0
  }
  
  if (sum(inds_hurts)>0) {
    res_hurts<-basic_lee_bound(leedata[inds_hurts,],...) 
    bounds_hurts<-GetBounds(res_hurts)
    s_hurts<-always_takers_share(leedata[inds_hurts,],...)
    
  } else {
    bounds_hurts<-c(0,0)
    s_hurts<-0
  }
  
  bounds<-(bounds_helps*s_helps+bounds_hurts*s_hurts)/(s_helps+s_hurts)
  
  return(list(lower_bound=bounds[1],upper_bound=bounds[2],
              bounds_helps=bounds_helps,
              bounds_hurts=bounds_hurts,
              s_helps=s_helps,
              s_hurts=s_hurts))
  
}
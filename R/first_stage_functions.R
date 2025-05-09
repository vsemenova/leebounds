### please explain what this function does
estimate_selection<-function(leedata,form=NULL,selection_function_name,variables_for_selection=NULL,names_to_include=c(),
                             treat_name="treat+",yname="selection",myweights=NULL,...) {
  ### read in data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  ### set variables for regression
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("treat","selection","outcome","X.Intercept.","(Intercept)","weights","prop1","prop0"))
  } else {
    variables_for_selection<-unique(setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)","weights","prop1","prop0")))
  }
  if (is.null(form)) {
    
    form<-as.formula(paste0("selection~(treat)*(", paste0(variables_for_selection,collapse="+"),")"))
  }
  
  print(form)
  
  if (is.null(myweights)) {
    myweights<-rep(1,dim(leedata)[1])
  }
  
  ## if rlassologit, select covariates
  if (selection_function_name=="rlassologit") {
    glm.fit<-rlassologit(form, leedata[,c("treat", "selection", variables_for_selection)],family="binomial",...)
    # select non-zero coefficients
    non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
    # names whose coefs are non-zero
    selected_names<-setdiff(names(non_zero_coefs),c("(Intercept)","treat"))
    # add manually selected features
    selected_names<-unique(c(selected_names,names_to_include))
    ## (optional): add raw variables after interactions
    #grep("treat:",selected_names,value=TRUE)
    
    # if treatment was dropped, make sure to re-run low-dim analysis with treatment
    ## this corresponds to assumption that Pr (S(1)=S(0)=1)<1 (not everyone is an always-taker)
    if (length(selected_names)>0) {
      form<-as.formula(paste0(yname,"~",treat_name,paste0(selected_names,collapse="+")))
    } else {
      form<-as.formula(paste0(yname,"~",treat_name))
    }
  }
  
  ### final stage is always logistic with low-dim covariates
  leedata$myweights<-myweights
  glm.postlasso<-glm( form,data=leedata[,c("treat", "selection", variables_for_selection)],family="binomial")
  
  return(glm.postlasso)
}

### please explain what this function does
predict_selection<-function(fit,leedata,...) {
  leedata_0treat<-leedata
  leedata_0treat$treat<-0
  
  leedata_1treat<-leedata
  leedata_1treat$treat<-1
  
  s.0.hat<-predict( fit,leedata_0treat,type="response")
  s.1.hat<-predict( fit,leedata_1treat,type="response")
  
  return(list(s.0.hat=s.0.hat,s.1.hat=s.1.hat))
}

### please explain what this function does
estimate_quantile<-function(leedata, form=NULL,
                            quantile_grid_size=0.01,
                            variables_for_outcome=NULL,
                            outcome_function_name="rq",...) {

  print(outcome_function_name)
  

  taus=seq(quantile_grid_size,1-quantile_grid_size,quantile_grid_size)
  estimated_quantiles_11<-matrix(0,dim(leedata)[1],length(taus))
  estimated_quantiles_10<-matrix(0,dim(leedata)[1],length(taus))
  
  outcome_data<-leedata[leedata$selection==1,]
  
  if (outcome_function_name=="rqLasso") {
    
    ## calculate data-driven lambda as in Belloni and Chernozhukov (2011)
    form_lasso<-as.formula(paste0("~treat*(", paste0(variables_for_outcome,collapse="+"),")"))
    ### covariate matrix with treatment interactions
    mydata<-model.matrix(form_lasso, outcome_data)
    ## covariates close to being degenerate
    singularity<-apply(mydata,2, var)
    count_sing <- sum(singularity <= 0.2)
    print("count of var that singularity < 0.1")
    print(count_sing)
    # select non-singular covariates
    # rq-lasso does not automatically drop degenerate covariates (need to do it instead)
    ## count of var that singularity <0.1
    ## examine identities var that can be dropped bc missing data/constant
    mydata<-mydata[,singularity>0.2]

    ##
    data_q<-as.data.frame(cbind(outcome=outcome_data$outcome,mydata))
    selected_covs_outcome<-c()  
    
    ## calculate Lambda as in Belloni and Chernozhukov uniformly over taus
    lambda<-lambda.BC(X=mydata,taus=taus)
    
    failures <- 0
    
    for (i in 1:(length(taus))) {
      tau<-taus[i]
      print("Estimating l1-regularized QR for tau=")
      print(tau)
    
  
      lambda_tau<- lambda * sqrt(tau*(1-tau))
      
      
      q_model<-quantreg::rq(outcome~0+.,data=data_q,tau=tau,method="lasso",lambda=lambda_tau)
      print("Selected coefficients")
      selected_coeffs <- names(q_model$coefficients[abs(q_model$coefficients) > 0.01])
      print(selected_coeffs)
      selected_covs_outcome<-c(selected_covs_outcome,names(q_model$coefficients[abs(q_model$coefficients)>0.01]) )

      leedata1<-leedata
      leedata1$treat <-1
      mydata1<-model.matrix(form_lasso, leedata1)
      data_q1<-as.data.frame(cbind(outcome=leedata$outcome,mydata1[,colnames(mydata)]))

      estimated_quantiles_11[,i]<-predict(q_model, data_q1)

      leedata0<-leedata
      leedata0$treat <-0
      mydata0<-model.matrix( form_lasso, leedata0)
      data_q0<-as.data.frame(cbind(outcome=leedata$outcome,mydata0[,colnames(mydata)]))

      estimated_quantiles_10[,i]<-predict(q_model, data_q0)
    }
    
    cat("Total failures:", failures, "\n")

    selected_covs_outcome=setdiff(unlist(sapply(selected_covs_outcome, function(str) strsplit(str,split="treat.") )), c("treat","(Intercept)"))
    print(selected_covs_outcome)
    
  } else {
    if (outcome_function_name=="rq") {
      
      for (i in 1:(length(taus))) {
        
        tau<-taus[i]
        print ("Estimating regular QR for tau=")
        print(tau)
        tau<-taus[i]
        
        data=outcome_data[,variables_for_outcome]
        
        singularity <- apply(data,2, var)
        data <- data[,singularity>0.001]
        variables_for_outcome <- names(data)
        
        ### set variables for regression
        variables_for_outcome<-unique(setdiff(variables_for_outcome,c("treat","selection")))
        
        if (is.null(form)) {
          form<-as.formula(paste0("outcome~(treat)*(", paste0(variables_for_outcome,collapse="+"),")"))
        }
        print(form)
        
        q_model<-quantreg::rq(form,data=outcome_data[,c("treat", "outcome",variables_for_outcome)],tau=tau)
        
        leedata1<-leedata
        leedata1$treat <-1
        
        estimated_quantiles_11[,i]<-predict(q_model, leedata1)
        
        leedata0<-leedata
        leedata0$treat <-0
        
        estimated_quantiles_10[,i]<-predict(q_model, leedata0)
      }
    }
      else {
        print ("No Quantile Regression method found")
      }
  }
  
  return(list(estimated_quantiles_11=estimated_quantiles_11,estimated_quantiles_10=estimated_quantiles_10))
}

impute_quantile_level<-function(quantile_table,p.0.hat,quantile_grid_size,min_wage, max_wage, ...) {
  
  taus=seq(quantile_grid_size,1-quantile_grid_size,quantile_grid_size)
  
  y.p0.hat<-rep(0,dim(quantile_table)[1])
  
  for (i in 1:length(taus)) {
    tau<-taus[i]
    inds<-abs(p.0.hat-tau)<=quantile_grid_size
    if(sum(inds)>0) {
      y.p0.hat[inds]<-quantile_table[inds,i]
    } 
  }
  
  inds<-(p.0.hat<=quantile_grid_size)
  if (sum(inds)>0) {
    y.p0.hat[inds]<-min_wage
  }
  
  inds<-(p.0.hat>=1-quantile_grid_size)
  if (sum(inds)>0) {
    y.p0.hat[inds]<-max_wage
  }
  
  return(y.p0.hat)
}

estimate_borderline_wage<-function(quantile_table,p.0.hat,min_wage=NULL,max_wage=NULL,...) {
  ### sort quantile table
  
  for (obs in 1:dim(quantile_table)[1]) {
    quantile_table[obs,]<-sort(quantile_table[obs,])
  }
  #
  y.p0.hat<-impute_quantile_level(quantile_table,p.0.hat,min_wage=min_wage, max_wage=max_wage, ...)
  y.1.p0.hat<-impute_quantile_level(quantile_table,1-p.0.hat,min_wage=min_wage, max_wage=max_wage,...)
  if (!is.null(min_wage)) {
    y.p0.hat<-sapply(y.p0.hat,max,min_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,max,min_wage)
  }
  if (!is.null(max_wage)) {
    y.p0.hat<-sapply(y.p0.hat,min,max_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,min,max_wage)
    
  }
  
  return(data.frame(y.p0.hat=y.p0.hat, y.1.p0.hat=y.1.p0.hat))
}


first_stage_wrapper<-function(leedata,
                              variables_for_outcome,
                              quantile_grid_size,
                              sort_quantiles=TRUE,
                              s.hat=NULL,
                              s_min=0.001,...) {
  sample_size<-dim(leedata)[1]
  weights<-leedata$weights
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  } 
  
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  p = dim(leedata)[2]-3
  
  if (!is.null(s.hat)) {
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
    
  } else {
    glm.fit<-estimate_selection(leedata,...)
    res_selection<-predict_selection(glm.fit,leedata,...)
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    p.0.star<-s.0.hat/s.1.hat
  }
  p.0.star<-(s.0.hat/s.1.hat)
  
  p.0.star<-round_p0(p.0.star,...)
  
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  y.hat<-data.frame(y.p0.hat=rep(NA,sample_size),y.1.p0.hat=rep(NA,sample_size))
  
  estimated_quantiles<-estimate_quantile(leedata=leedata,
                           # outcome_function_name=outcome_function_name,
                            variables_for_outcome=variables_for_outcome,
                            quantile_grid_size=quantile_grid_size,...)
  
  if (sum(inds_helps)>0) {
     
    y.hat.helps=estimate_borderline_wage(taus=taus,quantile_table=estimated_quantiles$estimated_quantiles_11[inds_helps,],p.0.hat=p.0.star[inds_helps],
                                         quantile_grid_size=quantile_grid_size,sort_quantiles=sort_quantiles,...)
    y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
    y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
  } else {
    estimated_quantiles_11<-NULL
  }
  
  if (sum(inds_hurts)>0) {
       y.hat.hurts=estimate_borderline_wage(taus=taus,quantile_table=estimated_quantiles$estimated_quantiles_10[inds_hurts,],p.0.hat=1/p.0.star[inds_hurts],
                                         quantile_grid_size=quantile_grid_size,...)
    y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
    y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
  }
  else {
    estimated_quantiles_10<-NULL
  }
  s.hat$s.0.hat<-sapply(s.hat$s.0.hat,max,s_min)
  s.hat$s.1.hat<-sapply(s.hat$s.1.hat,max,s_min)
  
  leedata$s.0.hat<-s.hat$s.0.hat
  leedata$s.1.hat<-s.hat$s.1.hat
  leedata$p.0.star<-p.0.star
  leedata$y.p0.hat<-y.hat$y.p0.hat
  leedata$y.1.p0.hat<-y.hat$y.1.p0.hat
  
  return(leedata)
}

### auxiliary functions for l1-regularized quantile regression
norm2n<- function(z){  sqrt(mean(z^2)) }



lambda.BC<- function(X, taus, R = 1000, c = 1.2, alpha = .1){
  
  set.seed(1)
  n <- nrow(X)
  sigs <- apply(X,2,norm2n)
  U <- matrix(runif(n * R),n)
  Lambda <- array(0, dim = c(ncol(X), R, length(taus)))
  
  for (i in 1:(length(taus))) {
    
    tau <- taus[i]
    R <- (t(X) %*% (tau - (U < tau)))/(sigs*sqrt(tau*(1-tau)))
    Lambda[, , i] <- abs(R)
    
  }
  ## this lambda.BC encodes equation (2.6) in https://arxiv.org/pdf/0904.2931.pdf (Belloni and Chernozhukov, 2011)
  r <- apply(Lambda, MARGIN = 2, FUN = max)
  
  lambda <- c * quantile(r, 1 - alpha) * sigs
  
  return(lambda)
}


# leebounds
Code associated with paper: "Generalized Lee Bounds" by Vira Semenova, Journal of Econometrics, forthcoming

This replication package contains files that implement sharp bounds on treatment effects in the presence of selection/nonresponse bias in randomized control trials. It includes basic Lee bounds  [Lee (2009)](https://academic.oup.com/restud/article-abstract/76/3/1071/1590707) and sharp bounds based on machine learning [Semenova(2023)](https://arxiv.org/abs/2008.12720). 

### Example of standard Lee Bounds

```
######### Compute basic Lee (2009) bounds for ATE in week 208 #########
leedata=data.frame(treat=JobCorps_baseline$TREATMNT.y,selection=JobCorps_employment$week_208,outcome=JobCorps_wages$week_208)
GetBounds(leebounds(leedata))
```

### Example of generalized Lee Bounds
```
######### Compute basic Lee (2009) bounds for ATE in week 208 #########
leedata=data.frame(treat=JobCorps_baseline$TREATMNT.y,selection=JobCorps_employment$week_208,outcome=JobCorps_wages$week_208)

### STEP1:  main arguments
Lee_data<-read.csv(paste0(my_path,"/JobCorps/Derived_Data/dataLee2009.csv"))

week<-90
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
leedata_cov<-prepare_leedata(week=week,Lee_data,covariate_data = Lee_data[,baseline_varnames])

### other parameters
N_rep=500
ci_alpha=0.05
quantile_grid_size=0.01
min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])

### trimming threshold
sample_size<-dim(leedata_cov)[1]
rhoN = 1*(sample_size)^{-1/4} * (log((sample_size)))^{-1}

leedata_cov$prop1=1/2*rep(1,sample_size)
leedata_cov$prop0=1/2*rep(1,sample_size)
p0=1-rhoN
p1=1+rhoN

res_wrapper<-ortho_leebounds(leedata=leedata_cov,
                           selection_function_name="glm",
                           variables_for_selection=baseline_varnames,
                           quantile_grid_size = quantile_grid_size,
                           variables_for_outcome=baseline_varnames,
                           min_wage=min_wage,
                           max_wage=max_wage,ortho=TRUE,
                           p0=p0, p1=p1,outcome_function_name="rq")


print("Estimated bounds: ortho_leebounds function")
print(GetBounds(res_wrapper))
print("Estimated always-takers' share")
print(res_wrapper$denom)
```

# Figures and Tables

```
module load R

cd JobCorps

cd STEP1_PREPARE_DATA

Rscript prepare_basic_data.R prepare_covaraites.R prepare_covaraites_short.R

```

To replicate Figures 1 and 2, run R files in STEP2_TEST_MONOTONICITY

```
module load R

cd JobCorps

cd STEP2_TEST_MONOTONICITY

Rscript Figures12.R test_monotonicity.R  print_Figures12.R
```

To replicate Tables 1 and 2, run R files in STEP3_ESTIMATE_BOUNDS

```
module load R

cd JobCorps

cd STEP3_ESTIMATE_BOUNDS

Rscript Tab_Lee.R Tab_Full.R Tab_Lasso.R Table1.R Table2.R
```


To replicate Figure 3, run R files in STEP4_FIGURE3

```
module load R

cd JobCorps

cd module load R

cd JobCorps

cd STEP4_FIGURE3

Rscript Figure3.R print_Figure3.R
```

# Support
Vira Semenova: semenovavira@gmail.com

# References
``Training, Wages, and Sample Selection: Estimating Sharp Bounds on Treatment Effects'' by David S. Lee, Review of Economic Studies (2009)76, 1071-1 102

``Generalized Lee Bounds'' by Vira Semenova, Journal of Econometrics, forthcoming

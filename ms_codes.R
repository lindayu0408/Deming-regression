library(dplyr)
library(ggplot2)
library(tidyr)
library(survival)
library(tableone)
library(ppcor)
library(correlation)
library(mcr)
library(gridExtra)

####################################################################################################

###############
# Table 1
###############

vars<-c("gender","mom_age_del","parity_ord","moeducat_ord","mom_race_eth","comorbidity_nodm","season_lmp2","year_3c","PCA","preg_pm25mass","preg_pm01mass")
factor_vars<-c("gender","parity_ord","moeducat_ord","mom_race_eth","comorbidity_nodm","year_3c")
t1<-CreateTableOne(data=data,vars=vars,factorVars=factor_vars)
print(t1)

t2<-CreateTableOne(data=data,vars=vars,factorVars=factor_vars,strata = "event_asd2")
print(t2)

####################################################################################################
data$preg_pm25mass_sd<-data$preg_pm25mass/sd(data$preg_pm25mass)
data$preg_pm01mass_sd<-data$preg_pm01mass/sd(data$preg_pm01mass)
data$preg_pm25_01mass<-data$preg_pm25mass-data$preg_pm01mass
data$preg_pm25_01mass_sd<-data$preg_pm25_01mass/sd(data$preg_pm25_01mass)

###################
# Table 2
##################

lm_fit_PM25_sd<-lm(preg_pm25mass_sd~preg_pm01mass_sd,data=data)
summary(lm_fit_PM25_sd) 
confint(lm_fit_PM25_sd)

data$fitted_pm25<-predict(lm_fit_PM25_sd)
data$resid_pm25<-resid(lm_fit_PM25_sd)

lm_fit_PM01_sd<-lm(preg_pm01mass_sd~preg_pm25mass_sd,data=data)
summary(lm_fit_PM01_sd) 
confint(lm_fit_PM01_sd)

data$fitted_pm01<-predict(lm_fit_PM01_sd)
data$resid_pm01<-resid(lm_fit_PM01_sd)

PM25_sd_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ preg_pm25mass_sd +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
summary(PM25_sd_fit1)

PM01_sd_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ preg_pm01mass_sd +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
summary(PM01_sd_fit1)

PM25_sd_decomp_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ fitted_pm25+resid_pm25 +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
summary(PM25_sd_decomp_fit1)

PM01_sd_decomp_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ fitted_pm01+resid_pm01 +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
summary(PM01_sd_decomp_fit1)

######################################################################################################

##################
# Figure 1
#################
lambda_l<-seq(0.5,2,0.1)

deming_res0<-rep(-9,2)

HR_res_xhat0<-rep(-9,3)
HR_res_xresid0<-rep(-9,3)
HR_res_yhat0<-rep(-9,3)
HR_res_yresid0<-rep(-9,3)

coef_res_xhat0<-rep(-9,6)
coef_res_xresid0<-rep(-9,6)
coef_res_yhat0<-rep(-9,6)
coef_res_yresid0<-rep(-9,6)

for (lambda in lambda_l){
  print(lambda)
  
  model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)
  printSummary(model1)
  deming_res0<-rbind(deming_res0,c(model1@para[1,1],model1@para[2,1]))
  
  resid<-getResiduals(model1)
  hat<-getFitted(model1)
  
  # raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
  x_hat<-hat$x_hat
  y_hat<-hat$y_hat
  x_resid<-resid$x
  y_resid<-resid$y
  optimal_resid<-resid$optimized
  
  preg_x_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ x_hat+x_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  summary(preg_x_fit1)
  
  fit_content<-summary(preg_x_fit1)
  xhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  xhat_coef<-fit_content$coefficients[1,]
  xresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  xresid_coef<-fit_content$coefficients[2,]
  
  HR_res_xhat0<-rbind(HR_res_xhat0,xhat_HRs)
  HR_res_xresid0<-rbind(HR_res_xresid0,xresid_HRs)
  coef_res_xhat0<-rbind(coef_res_xhat0,xhat_coef)
  coef_res_xresid0<-rbind(coef_res_xresid0,xresid_coef)
  
  preg_y_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ y_hat+y_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  summary(preg_y_fit1)
  
  fit_content<-summary(preg_y_fit1)
  yhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  yhat_coef<-fit_content$coefficients[1,]
  yresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  yresid_coef<-fit_content$coefficients[2,]
  
  HR_res_yhat0<-rbind(HR_res_yhat0,yhat_HRs)
  HR_res_yresid0<-rbind(HR_res_yresid0,yresid_HRs)
  coef_res_yhat0<-rbind(coef_res_yhat0,yhat_coef)
  coef_res_yresid0<-rbind(coef_res_yresid0,yresid_coef)
}

xhat_HRs_scale_df0<-data.frame(HR_res_xhat0[-1,])
xhat_coef_scale_df0<-data.frame(coef_res_xhat0[-1,])
xresid_HRs_scale_df0<-data.frame(HR_res_xresid0[-1,])
xresid_coef_scale_df0<-data.frame(coef_res_xresid0[-1,])
yhat_HRs_scale_df0<-data.frame(HR_res_yhat0[-1,])
yhat_coef_scale_df0<-data.frame(coef_res_yhat0[-1,])
yresid_HRs_scale_df0<-data.frame(HR_res_yresid0[-1,])
yresid_coef_scale_df0<-data.frame(coef_res_yresid0[-1,])
deming_res_scale_df0<-data.frame(deming_res0[-1,])

xhat_HRs_scale_df0$lambda<-lambda_l
xhat_coef_scale_df0$lambda<-lambda_l
xresid_HRs_scale_df0$lambda<-lambda_l
xresid_coef_scale_df0$lambda<-lambda_l
yhat_HRs_scale_df0$lambda<-lambda_l
yhat_coef_scale_df0$lambda<-lambda_l
yresid_HRs_scale_df0$lambda<-lambda_l
yresid_coef_scale_df0$lambda<-lambda_l
deming_res_scale_df0$lambda<-lambda_l

p1<-ggplot(xhat_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of error variances")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2)+ggtitle(expression("(a) Linear-transformed PM"[2.5]))

p2<-ggplot(xresid_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of error variances")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(b) Residuals of PM"[2.5]))

p3<-ggplot(yhat_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of error variances")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(c) Linear-transformed PM"[0.1]))

p4<-ggplot(yresid_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of error variances")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(d) Residuals of PM"[0.1]))

grid.arrange(p1, p2, p3, p4, ncol = 2)
####################################################################################################

###################
# eTable 1
###################

cor_df<-data%>%dplyr::select(preg_pm25mass_sd,preg_pm01mass_sd,preg_pm25_01mass_sd,fitted_pm25,fitted_pm01,resid_pm25,resid_pm01)
cor(cor_df)

##################
# eTable 2
##################
lambda<-2

model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
               mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)
printSummary(model1)

resid<-getResiduals(model1)
hat<-getFitted(model1)

# raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
x_hat<-hat$x_hat
y_hat<-hat$y_hat
x_resid<-resid$x
y_resid<-resid$y
optimal_resid<-resid$optimized

cor_df2<-data.frame(preg_pm25mass_sd=data$preg_pm25mass_sd,preg_pm01mass_sd=data$preg_pm01mass_sd,x_hat,y_hat,x_resid,y_resid)
cor(cor_df2)

####################
# eTable 3
###################

# deming regression
lambda_l<-seq(0.5,2,0.1)
intercept_l<-rep(-9,3)
slope_l<-rep(-9,3)

HR_res_xhat0<-rep(-9,3)
HR_res_xresid0<-rep(-9,3)
HR_res_yhat0<-rep(-9,3)
HR_res_yresid0<-rep(-9,3)

coef_res_xhat0<-rep(-9,6)
coef_res_xresid0<-rep(-9,6)
coef_res_yhat0<-rep(-9,6)
coef_res_yresid0<-rep(-9,6)

for (lambda in lambda_l){
  print(lambda)
  
  model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)
  printSummary(model1)
  intercept_l<-rbind(intercept_l,model1@para[1,])
  slope_l<-rbind(slope_l,model1@para[2,])
}
intercept_df<-data.frame(intercept_l[-1,])
slope_df<-data.frame(slope_l[-1,])

intercept_df$lambda<-lambda_l
slope_df$lambda<-lambda_l


# boot for WDeming regression

data2<-data%>%dplyr::select(preg_pm01mass_sd,preg_pm25mass_sd)
#lambda_l<-seq(1.1,2.0,0.1)
for (lambda in lambda_l){
  print(lambda)
  
  boot_p <- function(...){
    data<-data2
    R<-250
    ncpus<-4
    lambda<-lambda
    
    boot_fun <- function(data, indices) {
      # curVal <- get("counter", envir = env) + ncpus
      # setTxtProgressBar(get("progbar", envir = env), curVal)
      # assign("counter", curVal+1, envir = env)
      d <- data[indices,]
      model<- mcreg(d$preg_pm01mass_sd,d$preg_pm25mass_sd,error.ratio=lambda,method.reg="WDeming", method.ci="bootstrap",
                    mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE,nsamples = 1)
      printSummary(model)
      res<-c(model@glob.coef)
      return(res)
    }
    
    results <- boot(data=data, statistic=boot_fun, R=R)
    return(results)
  }
  
  start_time <- Sys.time()
  cl <- makeCluster(mc <- getOption("cl.cores",4))
  
  clusterSetRNGStream(cl, 1234)
  clusterEvalQ(cl, {library(dplyr)
    library(mcr)
    library(boot)
    library(parallel)})
  clusterExport(cl, varlist = c("data2","lambda"))
  
  ## to make this reproducible
  # clusterSetRNGStream(cl, 123)
  res <- do.call(c, parLapply(cl, seq_len(mc), boot_p))
  end_time <- Sys.time()
  stopCluster(cl)
  print(end_time - start_time)
  
  
  lambda_str<-gsub("\\.", "_", format(lambda, nsmall = 1))
  save(res,file=paste("Z:/APARKPSC/Xin/PM01/deming_res/wdeming_res/Wdeming_PM25_PM01_lambda",lambda_str,".Rdata",sep=""))
}

wdeming_intercept_l<-c(-9,-9,-9)
wdeming_slope_l<-c(-9,-9,-9)
for (lambda in lambda_l){
  lambda_str<-gsub("\\.", "_", format(lambda, nsmall = 1))
  load(paste("Z:/APARKPSC/Xin/PM01/deming_res/wdeming_res/Wdeming_PM25_PM01_lambda",lambda_str,".Rdata",sep=""))
  res$t0
  intercept_ci<-boot.ci(res,index=1, type = "perc")
  wdeming_intercept_l<-rbind(wdeming_intercept_l,c(res$t0[1],intercept_ci$percent[c(4,5)]))
  
  slope_ci<-boot.ci(res,index=2, type = "perc")
  wdeming_slope_l<-rbind(wdeming_slope_l,c(res$t0[2],slope_ci$percent[c(4,5)]))
}


####################
# eTable 4
###################
# Deming regression
lambda_l<-seq(0.5,2,0.1)

rsquared_l<-c(-9,-9)
for (lambda in lambda_l){
  print(lambda)
  
  model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)
  printSummary(model1)
  
  resid<-getResiduals(model1)
  hat<-getFitted(model1)
  
  # raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
  x_hat<-hat$x_hat
  y_hat<-hat$y_hat
  x_resid<-resid$x
  y_resid<-resid$y
  optimal_resid<-resid$optimized
  
  PM01_fit<-lm(data$preg_pm01mass_sd~x_hat)
  PM01_fit_content<-summary(PM01_fit)
  
  PM25_fit<-lm(data$preg_pm25mass_sd~y_hat)
  PM25_fit_content<-summary(PM25_fit)
  rsquared_l<-rbind(rsquared_l,c(PM25_fit_content$r.squared,PM01_fit_content$r.squared))
}

rsquared_df<-data.frame(rsquared_l[-1,])
colnames(rsquared_df)<-c("PM25R2","PM01R2")
rsquared_df$lambda<-lambda_l

# WDeming regression
lambda_l<-seq(0.5,2,0.1)

print(lambda)

model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="WDeming", method.ci="bootstrap",
               mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE,nsamples=1,iter.max = 100)
printSummary(model1)

model0<-mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
              mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)

hat<-getFitted(model0)
x_hat<-hat$x_hat
y_hat<-hat$y_hat

rsquared_l<-c(-9,-9)
for (lambda in lambda_l){
  print(lambda)
  
  model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="WDeming", method.ci="bootstrap",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE,nsamples=1)
  printSummary(model1)
  
  resid<-getResiduals(model1)
  hat<-getFitted(model1)
  
  # raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
  x_hat<-hat$x_hat
  y_hat<-hat$y_hat
  x_resid<-resid$x
  y_resid<-resid$y
  optimal_resid<-resid$optimized
  
  PM01_fit<-lm(data$preg_pm01mass_sd~x_hat)
  PM01_fit_content<-summary(PM01_fit)
  
  PM25_fit<-lm(data$preg_pm25mass_sd~y_hat)
  PM25_fit_content<-summary(PM25_fit)
  rsquared_l<-rbind(rsquared_l,c(PM25_fit_content$r.squared,PM01_fit_content$r.squared))
}

rsquared_df<-data.frame(rsquared_l[-1,])
colnames(rsquared_df)<-c("PM25R2","PM01R2")
rsquared_df$lambda<-lambda_l


#########################################################################################################
#####################
# eFigure 3
####################


lambda_l<-seq(0.5,2,0.1)

deming_res0<-rep(-9,2)

HR_res_xhat0<-rep(-9,3)
HR_res_xresid0<-rep(-9,3)
HR_res_yhat0<-rep(-9,3)
HR_res_yresid0<-rep(-9,3)

coef_res_xhat0<-rep(-9,6)
coef_res_xresid0<-rep(-9,6)
coef_res_yhat0<-rep(-9,6)
coef_res_yresid0<-rep(-9,6)

for (lambda in lambda_l){
  print(lambda)
  
  model1<- mcreg(data$preg_pm01mass_sd,data$preg_pm25_01mass_sd,error.ratio=lambda,method.reg="Deming", method.ci="analytical",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE)
  printSummary(model1)
  deming_res0<-rbind(deming_res0,c(model1@para[1,1],model1@para[2,1]))
  
  resid<-getResiduals(model1)
  hat<-getFitted(model1)
  
  # raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
  x_hat<-hat$x_hat
  y_hat<-hat$y_hat
  x_resid<-resid$x
  y_resid<-resid$y
  optimal_resid<-resid$optimized
  
  preg_x_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ x_hat+x_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  print(summary(preg_x_fit1))
  
  fit_content<-summary(preg_x_fit1)
  xhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  xhat_coef<-fit_content$coefficients[1,]
  xresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  xresid_coef<-fit_content$coefficients[2,]
  
  HR_res_xhat0<-rbind(HR_res_xhat0,xhat_HRs)
  HR_res_xresid0<-rbind(HR_res_xresid0,xresid_HRs)
  coef_res_xhat0<-rbind(coef_res_xhat0,xhat_coef)
  coef_res_xresid0<-rbind(coef_res_xresid0,xresid_coef)
  
  preg_y_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ y_hat+y_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  print(summary(preg_y_fit1))
  
  fit_content<-summary(preg_y_fit1)
  yhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  yhat_coef<-fit_content$coefficients[1,]
  yresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  yresid_coef<-fit_content$coefficients[2,]
  
  HR_res_yhat0<-rbind(HR_res_yhat0,yhat_HRs)
  HR_res_yresid0<-rbind(HR_res_yresid0,yresid_HRs)
  coef_res_yhat0<-rbind(coef_res_yhat0,yhat_coef)
  coef_res_yresid0<-rbind(coef_res_yresid0,yresid_coef)
}

xhat_HRs_scale_df0<-data.frame(HR_res_xhat0[-1,])
xhat_coef_scale_df0<-data.frame(coef_res_xhat0[-1,])
xresid_HRs_scale_df0<-data.frame(HR_res_xresid0[-1,])
xresid_coef_scale_df0<-data.frame(coef_res_xresid0[-1,])
yhat_HRs_scale_df0<-data.frame(HR_res_yhat0[-1,])
yhat_coef_scale_df0<-data.frame(coef_res_yhat0[-1,])
yresid_HRs_scale_df0<-data.frame(HR_res_yresid0[-1,])
yresid_coef_scale_df0<-data.frame(coef_res_yresid0[-1,])
deming_res_scale_df0<-data.frame(deming_res0[-1,])

xhat_HRs_scale_df0$lambda<-lambda_l
xhat_coef_scale_df0$lambda<-lambda_l
xresid_HRs_scale_df0$lambda<-lambda_l
xresid_coef_scale_df0$lambda<-lambda_l
yhat_HRs_scale_df0$lambda<-lambda_l
yhat_coef_scale_df0$lambda<-lambda_l
yresid_HRs_scale_df0$lambda<-lambda_l
yresid_coef_scale_df0$lambda<-lambda_l
deming_res_scale_df0$lambda<-lambda_l

p1<-ggplot(xhat_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(a) Linear-transformed PM"[0.1-2.5]))

p2<-ggplot(xresid_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(b) Residuals of PM"[0.1-2.5]))

p3<-ggplot(yhat_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(c) Linear-transformed PM"[0.1]))

p4<-ggplot(yresid_HRs_scale_df0,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(d) Residuals of PM"[0.1]))

grid.arrange(p1, p2, p3, p4, ncol = 2)
###################
# eFigure 4
###################
lambda_l<-seq(0.5,2,0.1)

deming_res<-rep(-9,2)

HR_res_xhat<-rep(-9,3)
HR_res_xresid<-rep(-9,3)
HR_res_yhat<-rep(-9,3)
HR_res_yresid<-rep(-9,3)

coef_res_xhat<-rep(-9,6)
coef_res_xresid<-rep(-9,6)
coef_res_yhat<-rep(-9,6)
coef_res_yresid<-rep(-9,6)

for (lambda in lambda_l){
  print(lambda)
  
  model2<- mcreg(data$preg_pm01mass_sd,data$preg_pm25mass_sd,error.ratio=lambda,method.reg="WDeming", method.ci="bootstrap",
                 method.bootstrap.ci="quantile",
                 mref.name = "scaled PM0.1mass", mtest.name = "scaled PM2.5mass", na.rm=TRUE,nsamples=1,iter.max = 10000)
  printSummary(model2)
  deming_res<-rbind(deming_res,c(model2@glob.coef))
  
  resid2<-getResiduals(model2)
  hat2<-getFitted(model2)
  
  # raw_resid<-data$preg_pm25mass-(deming_fit1[1]+deming_fit1[2]*data$preg_pm01mass)
  x_hat<-hat2$x_hat
  y_hat<-hat2$y_hat
  x_resid<-resid2$x
  y_resid<-resid2$y
  optimal_resid<-resid2$optimized
  
  preg_x_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ x_hat+x_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  summary(preg_x_fit1)
  
  fit_content<-summary(preg_x_fit1)
  xhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  xhat_coef<-fit_content$coefficients[1,]
  xresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  xresid_coef<-fit_content$coefficients[2,]
  
  HR_res_xhat<-rbind(HR_res_xhat,xhat_HRs)
  HR_res_xresid<-rbind(HR_res_xresid,xresid_HRs)
  coef_res_xhat<-rbind(coef_res_xhat,xhat_coef)
  coef_res_xresid<-rbind(coef_res_xresid,xresid_coef)
  
  preg_y_fit1<-coxph(Surv(fu_year_asd_age1_2, event_asd2) ~ y_hat+y_resid +pspline(year_2001,4) + mom_age_del + parity_ord + mom_race_eth + moeducat_ord +comorbidity_nodm +gender + PCA+season_lmp2,data=data)
  summary(preg_y_fit1)
  
  fit_content<-summary(preg_y_fit1)
  yhat_HRs<-fit_content$conf.int[1,c(1,3,4)]
  yhat_coef<-fit_content$coefficients[1,]
  yresid_HRs<-fit_content$conf.int[2,c(1,3,4)]
  yresid_coef<-fit_content$coefficients[2,]
  
  HR_res_yhat<-rbind(HR_res_yhat,yhat_HRs)
  HR_res_yresid<-rbind(HR_res_yresid,yresid_HRs)
  coef_res_yhat<-rbind(coef_res_yhat,yhat_coef)
  coef_res_yresid<-rbind(coef_res_yresid,yresid_coef)
}

xhat_HRs_scale_df<-data.frame(HR_res_xhat[-1,])
xhat_coef_scale_df<-data.frame(coef_res_xhat[-1,])
xresid_HRs_scale_df<-data.frame(HR_res_xresid[-1,])
xresid_coef_scale_df<-data.frame(coef_res_xresid[-1,])
yhat_HRs_scale_df<-data.frame(HR_res_yhat[-1,])
yhat_coef_scale_df<-data.frame(coef_res_yhat[-1,])
yresid_HRs_scale_df<-data.frame(HR_res_yresid[-1,])
yresid_coef_scale_df<-data.frame(coef_res_yresid[-1,])
deming_res_scale_df<-data.frame(deming_res[-1,])

xhat_HRs_scale_df$lambda<-lambda_l
xhat_coef_scale_df$lambda<-lambda_l
xresid_HRs_scale_df$lambda<-lambda_l
xresid_coef_scale_df$lambda<-lambda_l
yhat_HRs_scale_df$lambda<-lambda_l
yhat_coef_scale_df$lambda<-lambda_l
yresid_HRs_scale_df$lambda<-lambda_l
yresid_coef_scale_df$lambda<-lambda_l
deming_res_scale_df$lambda<-lambda_l

p1<-ggplot(xhat_HRs_scale_df,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(a) Linear-transformed PM"[2.5]))

p2<-ggplot(xresid_HRs_scale_df,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(b) Residuals of PM"[2.5]))

p3<-ggplot(yhat_HRs_scale_df,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(c) Linear-transformed PM"[0.1]))

p4<-ggplot(yresid_HRs_scale_df,aes(x=lambda,y=exp.coef.,ymin=lower..95,ymax=upper..95))+geom_point()+
  geom_errorbar(aes(ymin=lower..95, ymax=upper..95),width=0.05)+xlab("ratio of errors variance")+ylab("HR")+
  geom_hline(yintercept=1,linetype=2) +ggtitle(expression("(d) Residuals of PM"[0.1]))


grid.arrange(p1, p2, p3, p4, ncol = 2)

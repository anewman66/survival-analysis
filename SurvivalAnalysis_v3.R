#Install any of the following packages you don't have (not all are now necessary,
# but at some point they have been, so remove at your own risk)
BiocManager::install(c('epiDisplay',
                       'PoweR',
                       'devtools',
                       'stringi',
                       'survminer',
                       'dplyr',
                       'ggplot2','Rcpp',
                       'purrr',
                       'tidyr',
                       'lmtest',
                       'boot',
                       'colorspace',
                       'coxrt',
                       'pillar',
                       'devEMF',
                       'reshape'))

library('epiDisplay')
library('PoweR')
library('devtools')
library('stringi')
library('survminer')
library('dplyr')
library('ggplot2')
library('Rcpp')
library('purrr')
library('tidyr')
library('lmtest')
library('boot')
library('colorspace')
library('coxrt')
library('pillar')
library('devEMF')
library('reshape')
library('survival')

#Mac version:
setwd("/Users/alex/OneDrive - Teesside University/Student/2021/Newcastle")

survdata <- read.csv("Survival_Analysis_JB_v1.csv", header = TRUE, row.names=1)
surv_all <- survdata

#Subset for three main subgroups
#FAB/LMB Only
survdata <- surv_all
#survdata <- survdata[survdata$Include_Leuk_v2 == 'Yes',]

#Only the high risk group
#survdata <- survdata[!is.na(survdata$Check_TP53_High_Risk),]
str(survdata)

#Subset for only cases with data available during the pilot
survdata <- survdata[survdata$TP53.Abnormality. == 'No',]

#Check your data for any with only one level. Enter colnames below. 
drops <- c("TP53_MUT_C10orf71",
           "TP53_MUT_CCDC57",
           "TP53_MUT_SNCAIP",
           "TP53_MUT_TRPM4",
           "TP53_MUT_CFHR5",
           "TP53_MUT_HLTF",
           "TP53_MUT_HS3ST5")

survdata <- survdata[ , !(names(survdata) %in% drops)]

#sBL times (swap in OS.*.3yr or TTP.*.3yr depending on the analysis)
timeOS <- survdata$OS_3yr
eventOS <- survdata$OS_3yr_event

timeTTP <- survdata$TTP_3yr  
eventTTP <- survdata$TTP_3yr_event

timePFS <- survdata$PFS_3yr
eventPFS <- survdata$PFS_3yr_event

#Check the values are now in the survdata object
dim(survdata)
str(survdata) # This tells you what dimensions to put in the next command

#Subset the data to remove descriptor columns - removes an error when printing the univariate results.

ncol(survdata)
covariates <- survdata[,13:ncol(survdata)]
#covariates <- survdata[,19:ncol(survdata)]
covariates <- colnames(covariates)
covariates


univ_formulasOS <- sapply(covariates,
                          function(x) as.formula(paste('Surv(timeOS, eventOS)~', x)))

univ_formulasPFS <- sapply(covariates,
                           function(x) as.formula(paste('Surv(timePFS, eventPFS)~', x)))

univ_formulasTTP <- sapply(covariates,
                           function(x) as.formula(paste('Surv(timeTTP, eventTTP)~', x)))


#Running the Cox model on both OS and TTP datasets - NB: Ensure data isn't mright format
univ_modelsOS <- lapply(univ_formulasOS, function(x){coxph(x, data = survdata)})
univ_modelsPFS <- lapply(univ_formulasPFS, function(x){coxph(x, data = survdata)})
univ_modelsTTP <- lapply(univ_formulasTTP, function(x){coxph(x, data = survdata)})


# Extract data
univ_resultsOS <- lapply(univ_modelsOS,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           wald.test<-signif(x$wald["test"], digits=3)
                           beta<-signif(x$coef[1], digits=3);#coeficient beta
                           HR <-signif(x$coef[2], digits=3);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"], 3)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
#univ_OS_output <- data.frame(t(sapply(univ_resultsOS,c)))
#univ_OS_file <- do.call(rbind, univ_OS_output)

univ_resultsPFS <- lapply(univ_modelsPFS,
                          function(x){
                            x <- summary(x)
                            p.value<-signif(x$wald["pvalue"], digits=3)
                            wald.test<-signif(x$wald["test"], digits=3)
                            beta<-signif(x$coef[1], digits=3);#coeficient beta
                            HR <-signif(x$coef[2], digits=3);#exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"], 4)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                            res<-c(beta, HR, wald.test, p.value)
                            names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                          "p.value")
                            return(res)
                            #return(exp(cbind(coef(x),confint(x))))
                          })

#univ_PFS_output <- data.frame(t(sapply(univ_resultsPFS,c)))
#univ_PFS_file <- do.call(rbind, univ_PFS_output)

univ_resultsTTP <- lapply(univ_modelsTTP,
                          function(x){
                            x <- summary(x)
                            p.value<-signif(x$wald["pvalue"], digits=3)
                            wald.test<-signif(x$wald["test"], digits=3)
                            beta<-signif(x$coef[1], digits=3);#coeficient beta
                            HR <-signif(x$coef[2], digits=3);#exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"], 3)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                            res<-c(beta, HR, wald.test, p.value)
                            names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                          "p.value")
                            return(res)
                            #return(exp(cbind(coef(x),confint(x))))
                          })

#univ_TTP_output <- data.frame(t(sapply(univ_resultsTTP,c)))
#univ_TTP_file <- do.call(rbind, univ_TTP_output)

write.table(univ_OS_file, "univ_OS_output.txt",sep='\t',quote=F)
write.table(univ_PFS_file, "univ_PFS_output.txt",sep='\t',quote=F)
write.table(univ_TTP_end, "univ_TTP_output.txt",sep='\t',quote=F)





# Remember to rename output file here first!

resOS <- t(as.data.frame(univ_resultsOS, check.names = FALSE))
resOS <- as.data.frame(resOS)
resOS

write.table(resOS, "OS_corona1_stage.txt",sep='\t',quote = F)

resPFS <- t(as.data.frame(univ_resultsPFS, check.names = FALSE))
resPFS <- as.data.frame(resPFS)
resPFS

write.table(resPFS, "PFS_corona2_stage.txt",sep='\t',quote = F)

resTTP <- t(as.data.frame(univ_resultsTTP, check.names = FALSE))
resTTP <- as.data.frame(resTTP)
resTTP

write.table(resTTP, "TTP_corona1_stage.txt",sep='\t',quote = F)



#Cohort survival time - OS
surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ Any_TP53, data = survdata)
resOS <- summary(fit, times = c(3))
save.df.OS <- as.data.frame(resOS[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
save.df.OS

#Cohort survival time - TTP
surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ 1, data = survdata)
resTTP <- summary(fit, times = c(3))
save.df.TTP <- as.data.frame(resTTP[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
save.df.TTP

#Cohort survival time - PFS
surv <- Surv(timePFS,eventPFS)
fit <- survfit(surv ~ Any_TP53, data = survdata)
resPFS <- summary(fit, times = c(3))
save.df.PFS <- as.data.frame(resPFS[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
save.df.PFS

#Cohort survival time - PFS
surv <- Surv(timePFS,eventPFS)
fit <- survfit(surv ~ Any_TP53, data = survdata)
resPFS <- summary(fit, times = c(3))

surv <- Surv(timePFS,eventPFS)
fit <- survfit(surv ~ Check_TP53_High_Risk, data = survdata)
resPFS <- summary(fit, times = c(3))

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ Check_TP53_High_Risk, data = survdata)
resOS <- summary(fit, times = c(3))


#Median OS calulation when survival doesn't go below 50%.
surv <- Surv(timeOS,eventOS)
fitOS <- survfit(surv ~ survdata$OS_3yr_event ==1 , data = survdata)
fitOS
plot(fitOS)


#Median TTP calulation when survival doesn't go below 50%.
surv <- Surv(timePFS,eventPFS)
fitTTP <- survfit(surv ~ survdata$PFS_3yr_event ==1 , data = survdata)
fitTTP
plot(fitTTP)


# Multivariate for the RR - forward selection.
TP53_bi_fit <- coxph(Surv(timeTTP, eventTTP) ~
                       survdata$Biallelic_vs_normal +
                       survdata$CNS +
                       survdata$BM 
                     ,data = survdata)
summary(TP53_bi_fit)

TP53_mono_fit <- coxph(Surv(timeTTP, eventTTP) ~
                         survdata$Mono_vs_normal +
                         survdata$CNS +
                         survdata$BM 
                       ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timeTTP, eventTTP) ~
                        survdata$Any_TP53 +
                        survdata$CNS +
                        survdata$BM 
                      ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timeTTP, eventTTP) ~
                         survdata$Mono_vs_other +
                         survdata$Biallelic_vs_other +
                         survdata$CNS +
                         survdata$BM 
                       ,data = survdata)
summary(TP53_both_fit)

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$TTP.Event.3yr ==1 , data = survdata)
fitTTP
plot(fitTTP)

#See if Bone Marrow, CSF, CNS, Stage etc add to model
sFitClin1 <- coxph(Surv(time, event) ~
                     survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                     survdata$X17p.CNN.LOH +
                     survdata$X17q.CNN.LOH..manual...10mb. +
                     survdata$Bone.marrow.involvement +
                     survdata$CNS.involvement +
                     survdata$Disease.Stage..St.Jude.
                   ,data = survdata)
summary(sFitClin1)
lrtest(sFitClin1,sFit10)

lrtest(sFit2,sFit1)
lrtest(sFit3,sFit1)
lrtest(sFit4,sFit1)
lrtest(sFit5,sFit1)
lrtest(sFit6,sFit2)

#Multivariate for the RR - backward selection.

bFit1 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit1)
#Removing one factor each time
bFit2 <- coxph(Surv(time, event) ~
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit2)
lrtest(bFit2,bFit1)

bFit3 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit3)
lrtest(bFit3,bFit1)

bFit4 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit4)
lrtest(bFit4,bFit1)

bFit5 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit5)
lrtest(bFit5,bFit1)

bFit6 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain + 
                 survdata$X17p.CNN.LOH + 
                 survdata$X17q.CNN.LOH..manual...10mb. + 
                 survdata$Age..if.known. 
               ,data = survdata)
summary(bFit6)
lrtest(bFit6,bFit1)


summary(sFit6)
HR <- round(exp(coef(sFit)), 2)
CI <- round(exp(confint(sFit)), 2)
P <- round(coef(summary(sFit))[,5], 3)
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")
# Bind columns together as dataset
table2 <- as.data.frame(cbind(HR, CI, P))
table2
#  #   #write.csv(table2, "RR_MultivariateForTier3_analysis2.csv")

######## BIALLELIC STUFF ########
TP53Fit0 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain
                  ,data = survdata)
summary(TP53Fit0)

TP53Fit1 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH
                  ,data = survdata)
summary(TP53Fit1)
lrtest(TP53Fit1,TP53Fit0)

TP53Fit2 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17q.CNN.LOH..manual...10mb.
                  ,data = survdata)
summary(TP53Fit2)
lrtest(TP53Fit2,TP53Fit0)

TP53Fit3 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$TP53.Biallelic
                  ,data = survdata)
summary(TP53Fit3)
lrtest(TP53Fit3,TP53Fit0)

TP53Fit4 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$TP53.Mutation..Y.N.
                  ,data = survdata)
summary(TP53Fit4)
lrtest(TP53Fit4,TP53Fit0)

TP53Fit5 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$X17q.CNN.LOH..manual...10mb.
                  ,data = survdata)
summary(TP53Fit5)
lrtest(TP53Fit5,TP53Fit1)

TP53Fit6 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$TP53.Biallelic
                  ,data = survdata)
summary(TP53Fit6)
lrtest(TP53Fit6,TP53Fit1)

TP53Fit7 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$TP53.Mutation..Y.N.
                  ,data = survdata)
summary(TP53Fit7)
lrtest(TP53Fit7,TP53Fit1)


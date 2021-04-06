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

### Part 1: Setup working directory and files. ###

#Remember if you paste the path in from Windows Explorer to switch the slashes from \ to /

#Specify working directory.
setwd("/Users/alex/OneDrive - Teesside University/Survival")

survdata <- read.csv("Form_eBL_Data.csv", header = TRUE, row.names=1)
surv_all <- survdata

View(survdata)

survdata <- surv_all

time <- survdata$RR.time..mo.
event <- survdata$RR.event

#Check the values are now in the survdata object
dim(survdata)
str(survdata) # This tells you what dimensions to put in the next command

#Subset the data to remove descriptor columns - removes an error when printing the univariate results.
ncol(survdata)
covariates <- survdata[,8:ncol(survdata)] #Note - change the number here to say where your data actually starts.
covariates <- colnames(covariates)
covariates

### Part 2: Run univariate analysis ###

#The following steps may provide warnings, but shouldn't provide errors.

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, event)~', x)))

#Running the Cox model on both OS and TTP datasets - NB: Ensure data isn't mright format.
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = survdata)})

# Extract data - if you have non-dichotomised input data (i.e. more than 2 possible outcomes) then this won't be converted to file properly. 
# Print "univ_results" to see this and just copy to a spreadsheet from the console.                     
univ_results <- lapply(univ_models,
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

# Remember to rename output file here first!

res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
res

write.table(res, "results.txt",sep='\t',quote = F)



                        
                        
                        
                        
### Part 3: Survival Metrics - Cohort Survival time, etc. ###

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


surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ Check_TP53_High_Risk, data = survdata)
resOS <- summary(fit, times = c(3))

surv <- Surv(timePFS,eventPFS)
fit <- survfit(surv ~ InterB_Group, data = survdata)
resPFS <- summary(fit, times = c(3))

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ InterB_Group, data = survdata)
resOS <- summary(fit, times = c(3))








fit <- survfit(surv ~ Int_High_Other, data = survdata)
resPFS <- summary(fit, times = c(3))

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ Any_TP53, data = survdata)
resPFS <- summary(fit, times = c(3))

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

### Part 4: Multivariate Analysis.
                           
                           
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

# Multivariate for the OS - forward selection.
TP53_bi_fit <- coxph(Surv(timeOS, eventOS) ~
                       survdata$Biallelic_vs_normal +
                       survdata$CNS +
                       survdata$BM 
                     ,data = survdata)
summary(TP53_bi_fit)

TP53_mono_fit <- coxph(Surv(timeOS, eventOS) ~
                         survdata$Mono_vs_normal +
                         survdata$CNS +
                         survdata$BM 
                       ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timeOS, eventOS) ~
                        survdata$Any_TP53 +
                        survdata$CNS +
                        survdata$BM +
                        survdata$LDH_high_v2 +
                        survdata$Stage_split
                      ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timeOS, eventOS) ~
                         survdata$Mono_vs_other +
                         survdata$Biallelic_vs_other +
                         survdata$CNS +
                         survdata$BM +
                         survdata$LDH_high_v2 +
                         survdata$Stage_split
                       ,data = survdata)
summary(TP53_both_fit)



# Multivariate for the RR - forward selection.
TP53_bi_fit <- coxph(Surv(timePFS, eventPFS) ~
                       survdata$Biallelic_vs_normal +
                       survdata$CNS +
                       survdata$BM + 
                       survdata$LDH_high
                     ,data = survdata)
summary(TP53_bi_fit)

TP53_LDH <- coxph(Surv(timePFS, eventPFS) ~
                       survdata$Any_TP53 +
                       survdata$LDH_high_v2
                     ,data = survdata)
summary(TP53_LDH)

TP53_mono_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_normal +
                         survdata$CNS +
                         survdata$BM + 
                         survdata$LDH_high
                       ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$Any_TP53 +
                        survdata$CNS +
                        survdata$BM +
                        survdata$LDH_high_v2 +
                        survdata$Stage_split
                      ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_other +
                         survdata$Biallelic_vs_other +
                         survdata$CNS +
                         survdata$BM +
                         survdata$LDH_high_v2 +
                         survdata$Stage_split
                       ,data = survdata)
summary(TP53_both_fit)

# NO LDH - PFS 
TP53_bi_fit <- coxph(Surv(timePFS, eventPFS) ~
                       survdata$Biallelic_vs_normal +
                       survdata$CNS +
                       survdata$BM +
                       survdata$LDH_high + 
                       survdata$Stage_split
                     ,data = survdata)
summary(TP53_bi_fit)

TP53_mono_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_normal +
                         survdata$CNS +
                         survdata$BM +
                         survdata$LDH_high + 
                         survdata$Stage_split
                       ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$Any_TP53 +
                        survdata$BM +
                        survdata$CNS +
                        survdata$LDH_high +
                        survdata$Treatment_BvC +
                        survdata$Stage_split
                      ,data = survdata)
summary(TP53_any_fit)

TP53_any_fit <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$Any_TP53 +
                        survdata$LDH_high +
                        survdata$Stage_split
                      ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_other +
                         survdata$Biallelic_vs_other +
                         survdata$LDH_high +
                         survdata$Stage_split
                       ,data = survdata)
summary(TP53_both_fit)

# With age - PFS

# NO LDH - PFS 
TP53_bi_fit <- coxph(Surv(timePFS, eventPFS) ~
                       survdata$Biallelic_vs_normal +
                       survdata$CNS +
                       survdata$BM +
                       survdata$Age
                     ,data = survdata)
summary(TP53_bi_fit)

TP53_mono_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_normal +
                         survdata$CNS +
                         survdata$BM +
                         survdata$Age
                       ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$Any_TP53 +
                        survdata$CNS +
                        survdata$BM +
                        survdata$Age
                      ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timePFS, eventPFS) ~
                         survdata$Mono_vs_other +
                         survdata$Biallelic_vs_other +
                         survdata$CNS +
                         survdata$BM +
                         survdata$Age
                       ,data = survdata)
summary(TP53_both_fit)




#Complexity

fit_13q <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$X13q +
                        survdata$CNS +
                        survdata$BM +
                        survdata$LDH_high_v2 + 
                        survdata$Stage_split
                      ,data = survdata)
summary(fit_13q)

fit_11q <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$X11q +
                        survdata$CNS +
                        survdata$BM +
                        survdata$LDH_high_v2 + 
                        survdata$Stage_split
                      ,data = survdata)
summary(fit_11q)

fit_1q <- coxph(Surv(timePFS, eventPFS) ~
                        survdata$X1q +
                        survdata$CNS +
                        survdata$BM +
                        survdata$LDH_high_v2 + 
                        survdata$Stage_split
                      ,data = survdata)
summary(fit_1q)


#OS Complexity
fit_13q <- coxph(Surv(timeOS, eventOS) ~
                   survdata$X13q +
                   survdata$CNS +
                   survdata$BM +
                   survdata$LDH_high_v2 + 
                   survdata$Stage_split
                 ,data = survdata)
summary(fit_13q)

fit_11q <- coxph(Surv(timeOS, eventOS) ~
                   survdata$X11q +
                   survdata$CNS +
                   survdata$BM +
                   survdata$LDH_high_v2 + 
                   survdata$Stage_split
                 ,data = survdata)
summary(fit_11q)

fit_1q <- coxph(Surv(timeOS, eventOS) ~
                  survdata$X1q +
                  survdata$CNS +
                  survdata$BM +
                  survdata$LDH_high_v2 + 
                  survdata$Stage_split
                ,data = survdata)
summary(fit_1q)










fit_13q <- coxph(Surv(timeOS, eventOS) ~
                   survdata$X13q +
                   survdata$CNS +
                   survdata$LDH_high_v2
                   survdata$BM 
                 ,data = survdata)
summary(fit_13q)

fit_11q <- coxph(Surv(timeOS, eventOS) ~
                   survdata$X11q +
                   survdata$CNS +
                   survdata$LDH_high_v2
                   survdata$BM 
                 ,data = survdata)
summary(fit_11q)

fit_1q <- coxph(Surv(timeOS, eventOS) ~
                  survdata$X1q +
                  survdata$CNS +
                  survdata$LDH_high_v2
                  survdata$BM 
                ,data = survdata)
summary(fit_1q)

# eFit2 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$GenderF
#                ,data = survdata)
# summary(eFit2)
# lrtest(eFit2,eFit1)
# eFit3 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$MCL1.Gain_J19.
#                ,data = survdata)
# summary(eFit3)
# lrtest(eFit3,eFit1)
# eFit4 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$MCL1.Gain_J19. +
#                  survdata$GenderF
#                ,data = survdata)
# summary(eFit4)
# lrtest(eFit4,eFit3)
# 
# 
# sFit5 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$Age..if.known.
#                ,data = survdata)
# summary(sFit5)
# lrtest(sFit5,sFit1)
# sFit6 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17p.CNN.LOH +
#                  survdata$X17q.CNN.LOH..manual...10mb.
#                ,data = survdata)
# summary(sFit6)
# lrtest(sFit6,sFit3)
# sFit7 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17q.CNN.LOH..manual...10mb.+
#                  survdata$No_abn_50kb
#                ,data = survdata)
# summary(sFit7)
# lrtest(sFit7,sFit3)
# sFit8 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17q.CNN.LOH..manual...10mb. +
#                  survdata$Age..if.known.
#                ,data = survdata)
# summary(sFit8)
# lrtest(sFit8,sFit3):
#   
#   sFit9 <- coxph(Surv(time, event) ~
#                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                    survdata$X17p.CNN.LOH +
#                    survdata$X17q.CNN.LOH..manual...10mb. +
#                    survdata$No_abn_50kb
#                  ,data = survdata)
# summary(sFit9)
# lrtest(sFit9,sFit6)
# 
# sFit10 <- coxph(Surv(time, event) ~
#                   survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                   survdata$X17p.CNN.LOH +
#                   survdata$X17q.CNN.LOH..manual...10mb. +
#                   survdata$Age..if.known.
#                 ,data = survdata)
# summary(sFit10)
# lrtest(sFit10,sFit6)
# 
# sFit11 <- coxph(Surv(time, event) ~
#                   survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                   survdata$X17p.CNN.LOH +
#                   survdata$X17q.CNN.LOH..manual...10mb. +
#                   survdata$Age..if.known. +
#                   survdata$No_abn_50kb
#                 ,data = survdata)
# summary(sFit11)
# lrtest(sFit11,sFit10)

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

#Survival curves

#Cohort Plot
surv <- Surv(time,event)
fit <- survfit(surv ~ 1, data = survdata)
names(fit$strata)
summary(fit)
plot(fit)
pdf ("OS_Cohort_Malawi_all.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Cohort.",
           legend.labs = c("Cohort"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           #pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = T,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

##########################################################################################

#Any TP53 TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Any_TP53, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Any_TP53.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Abnormality",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#17p CNN-LOH TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$X17p_CNN.LOH, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_17p_CNN-LOH.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p CNN-LOH",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# 17p Loss

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$X17p_Loss, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_17p_Loss.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p Loss",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#TP53 Mutation TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Mutation, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_TP53_mutation.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Mutation",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Biallelic vs normal TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Biallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Biallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#Monoallelic vs others TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Monoallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Monoallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Monoallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# TP53 Split Plot TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic.Monoallelic.Normal_Exc_Silent, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_TP53_Split.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Status",
           legend.labs = c("Biallelic","Monoallelic","Normal"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = F,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","black"),
           linetype = c("solid","solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

############################################ OS ###

#Any TP53 OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Any_TP53, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Any_TP53.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Abnormality",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#17p CNN-LOH OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$X17p_CNN.LOH, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_17p_CNN-LOH.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p CNN-LOH",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# 17p Loss OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$X17p_Loss, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_17p_Loss.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p Loss",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#TP53 Mutation OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Mutation, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_TP53_mutation.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Mutation",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Biallelic vs normal OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Biallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Biallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#Monoallelic vs others OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Monoallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Monoallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Monoallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# TP53 Split Plot OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Biallelic.Monoallelic.Normal_Exc_Silent, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_TP53_Split.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Status",
           legend.labs = c("Biallelic","Monoallelic","Normal"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = F,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","black"),
           linetype = c("solid","solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()





#Cohort Plot
surv <- Surv(time,event)
fit <- survfit(surv ~ 1, data = survdata)
names(fit$strata)
summary(fit)
plot(fit)
pdf ("OS_Cohort_Malawi_all.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Cohort.",
           legend.labs = c("Cohort"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           #pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = T,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Plot Gender

surv <- Surv(time,event)
fit <- survfit(surv ~ survdata$Gender, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Gender.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Gender",
           legend.labs = c("Female","Male"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# InterB








#Log Rank Test
mir <-survdiff(Surv(timeOS,eventOS) ~ MIR17HG.Gain,data=survdata)
pvalue(mir)

survdata$Int_High_Low <- NA
survdata$Int_High_Low[survdata$InterB_Group == "High"] <- "High"
survdata$Int_High_Low[survdata$InterB_Group == "Low"] <- "Low"
HvL <- survdiff(Surv(timePFS,eventPFS) ~ Int_High_Low,data=survdata)
pvalue(HvL)

survdata$Int_High_Low <- NA
survdata$Int_High_Low[survdata$InterB_Group == "High"] <- "High"
survdata$Int_High_Low[survdata$InterB_Group == "Low"] <- "Low"
HvL <- survdiff(Surv(timeOS,eventOS) ~ Int_High_Low,data=survdata)
pvalue(HvL)

survdata$Int_High_Int[survdata$InterB_Group == "High"] <- "High"
survdata$Int_High_Int[survdata$InterB_Group == "Intermediate"] <- "Int"
HvI <- survdiff(Surv(timePFS,eventPFS) ~ Int_High_Int,data=survdata)
pvalue(HvI)

survdata$Int_High_Int[survdata$InterB_Group == "High"] <- "High"
survdata$Int_High_Int[survdata$InterB_Group == "Intermediate"] <- "Int"
HvI <- survdiff(Surv(timeOS,eventOS) ~ Int_High_Int,data=survdata)
pvalue(HvI)

survdata$Int_High_Other[survdata$InterB_Group == "High"] <- "High"
survdata$Int_High_Other[survdata$InterB_Group == "Intermediate"] <- "Other"
survdata$Int_High_Other[survdata$InterB_Group == "Low"] <- "Other"
HvI <- survdiff(Surv(timePFS,eventPFS) ~ Int_High_Int,data=survdata)
pvalue(HvI)

pvalue <- function(x, ...) UseMethod("pvalue")
pvalue.survdiff <- function (x, ...) 
{
  if (length(x$n) == 1) {
    df <- 1
    pval <- pchisq(x$chisq, 1, lower.tail = FALSE)
  } else {
    if (is.matrix(x$obs)) {
      otmp <- rowSums(x$obs)
      etmp <- rowSums(x$exp)
    } else {
      otmp <- x$obs
      etmp <- x$exp
    }
    df <- sum(etmp > 0) - 1
    pval <- pchisq(x$chisq, df, lower.tail = FALSE)
  }
  list(chisq = x$chisq, p.value = pval, df = df)
}


aba <- survdiff(Surv(timePFS,eventPFS) ~ Check_TP53_High_Risk,data=survdata)
aba <- survdiff(Surv(timeOS,eventOS) ~ Check_TP53_High_Risk,data=survdata)

# Section below is for the old favourite of dashed and solid black and white lines in KMs.
# surv <- Surv(time,event)
# fit <- survfit(surv ~ survdata$v9_GISTIC2_chr3.195718766.198022430..Gain, data = survdata)
# names(fit$strata)
# #names(fit$strata) <- c("No 17p","Yes 17p")
# plot(fit)
# pdf ("RR_3q29_gain_all.file")
# #ggsurvplot(fit,font.legend = c(6, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
# ggsurvplot(fit,
#            # Change legends: title & labels
#            legend.title = "TP53 Biallelic",
#            legend.labs = c("No", "Yes"),
#            xlab="Overall Survival in years",
#            # Add p-value and tervals
#            pval = FALSE,
#            
#            conf.int = FALSE,
#            # Add risk table
#            risk.table = TRUE,
#            tables.height = 0.2,
#            tables.theme = clean_theme(),
#            
#            # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
#            # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
#            palette = c("#0f0e0c", "#0f0e0c"),
#            linetype = c("solid","dashed"),
#            ggtheme = theme_survminer() # Change ggplot2 theme
# )
# dev.off()



surv <- Surv(time,event)
# fit <- survfit(surv ~ 1 , data = survdata)
# res <- summary(fit, times = c(1,2,3))
# save.df <- as.data.frame(res[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
# save.df
# write.csv(save.df, "RR_tier3_KmEstimates.csv")

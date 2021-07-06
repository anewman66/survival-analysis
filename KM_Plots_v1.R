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

#Remember if you paste the path in from Windows Explorer to switch the slashes from \ to /

setwd("")
full_input <- read.csv("input.csv", header = TRUE, row.names=1)

### Set your cohort

#All
survdata <- full_input
dim(survdata)

#sBL times (swap in OS.*.3yr or TTP.*.3yr depending on the analysis)
timeOS <- survdata$OS_3yr
eventOS <- survdata$OS_3yr_event

timePFS <- survdata$PFS_3yr
eventPFS <- survdata$PFS_3yr_event

timeTTP <- survdata$TTP_3yr  
eventTTP <- survdata$TTP_3yr_event

### Check the data is showing what you want ###

dim(survdata)
str(survdata)


### Set up what format you want for the KM plots

pval <- TRUE
conf.int <- F


#Run this whole block at once.
# Change "Gene" in "surv ~ survdata$Gene" to be the variable you want to plot.

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Gene, data = survdata)
names(fit$strata)


plot(fit)
pdf ("filename_output.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Gene of interest",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Time to Progression Survival",
           # Add p-value and tervals
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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


surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
fit

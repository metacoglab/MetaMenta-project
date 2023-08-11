## cleans up and analyses online prolific data (experiment 1)
## Original code by elisavanderplasATgmail.com
## Updated August 2023 by stephen.flemingATucl.ac.uk

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
require(ggplot2)
require(plyr)
require(readr)#to easily read csv's
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe

### 1) LOAD DATA

dat=NULL #initiate vars to load
fullIDs=NULL
sessions = c('v38', 'v39', 'v40', 'v41', 'v42', 'v43')
#load in all task data
for (i in 1:length(sessions)){ 
  currentData = paste('data_exp_12022-',sessions[i],sep="")
  dataDir = "~/Dropbox/InPrep/Autism/Code/MetaMenta-project-github/Data/Exp1/" 
  Dir = paste(dataDir,currentData, '/', sep="")
  files = c('_task-pf6t', '_task-yzt9')
  for (j in 1:length(files)){
    data = read.csv2(paste(Dir, currentData, files[j],'.csv',sep=""), header=T, sep=",", )
    dat = rbind(dat, data)
  }
  IDs = read.csv2(paste(Dir, 'total_IDs_', sessions[i], '.csv', sep = ""), header = F, sep = ",", )
  fullIDs = rbind(fullIDs, IDs)
}

#get the demographics
setwd(dataDir)
subData <- read_csv('data.csv')

## 2) EXTRACT DATA AND PUT IN BIG ARRAY
#prepare the meta-task data
fullIDs = t(fullIDs)
bigData = NULL
for (s in 1:length(fullIDs))
{
  subj_dat = dat[dat$Participant.Private.ID==fullIDs[s],] ##load variables for specific subject
  conftask_dat=subj_dat[subj_dat$Task_type=="simpleperceptual",]##select trials from the confidence sub-task
 
  vistrials=conftask_dat[conftask_dat$label=="responsePerceptual",]##select initial binary decision (left/right) trials from the confidence sub-task
  conftrials=conftask_dat[conftask_dat$label=="confidencerating",]##select subsequent confidence rating from the confidence sub-task

  logRT = scale(log(vistrials$Reactiontime))
  conf = round(as.numeric(conftrials$confidence_rating)*100)/100
  keypress=vistrials$key_press
  acc=vistrials$correct 
  acc[is.na(acc)]=0##accuracy==1: correct, accuracy==0: wrong
  acc = acc-0.5##accuracy==0.5: correct, accuracy==-0.5: wrong
  
  #have to compute objectively correct answer because js script doesn't give that 
  dir=rep(1, length(acc))
  for (t in 1:length(acc)){
    if (acc[t] ==1 & keypress[t] == 87){
      dir[t]= -1}
    else if (acc[t] == 0 && keypress[t]==69){
      dir[t] = -1}
    }##correct and chose left, dir == -1 (left) or wrong and chose right, dir == -1 (left)
    
  #get all vars behind each other per subject
  subj = rep(s, length(acc))
  AQ=rep(subData$AQ10[s], length(acc))
  RAADS = rep(subData$RAADS[s], length(acc))
  AQ_comm = rep(subData$AQ10_C[s], length(acc))
  RAADS_menta = rep(subData$RAADS_MENTA[s], length(acc))
  MCQ_cat = rep(subData$MCQ_cat[s], length(acc))
  MCQ_emos = rep(subData$MCQ_feelings[s], length(acc))
  subData1 = data.frame("subj"=subj, "dir"=dir,"AQ10"=AQ,"RAADS14"= RAADS,"AQ10_comm"=AQ_comm,"RAADS14_menta"= RAADS_menta, "acc"=acc, "conf"=conf, "logRT"=logRT, "MCQ_emo"=MCQ_emos, "MCQ_cat"=MCQ_cat)
    
  #add to larger file 
  bigData = rbind(bigData, subData1)
}
# Factors
bigData$subj <- factor(bigData$subj)
# skip rows with NaNs
bigData_clean <- na.omit(bigData)

## 3) COMPUTE INTERACTIONS BETWEEN RT-CONF EFFECT AND MCQ SCORES 

## 3a) MCQ-cat 

confModel_noMCQ = lmer(conf ~ acc + logRT + acc * logRT + (1 + acc + logRT|subj), data=bigData_clean
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

confModel_wMCQ = lmer(conf ~ MCQ_cat*(acc + logRT + acc * logRT) + (1 + acc + logRT|subj), data=bigData_clean
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wMCQ)
print(summary(confModel_wMCQ))
print(Anova(confModel_wMCQ, type = 3))
coef(summary(confModel_wMCQ)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wMCQ))) 

## check if including MCQ as interaction improved the fit of the model, line 257
anova(confModel_noMCQ,confModel_wMCQ) 

# get median MCQ
MCQ_cat_med <- median(bigData_clean$MCQ_cat)
bigData_clean$MCQcat_bi = cut(as.numeric(bigData_clean$MCQ_cat),breaks=c(0,MCQ_cat_med,max(bigData_clean$MCQ_cat)),include.lowest=T, labels=c("low", "high"))
bigData_clean$MCQcat_bi <- factor(bigData_clean$MCQcat_bi)

#make a nice figure of conf~RTxMCQ interaction 
ggplot(bigData_clean, aes(x=logRT, y=conf, colour=MCQcat_bi)) + 
  geom_count() + 
  scale_color_manual(values=c("salmon", "turquoise3")) +
  geom_point(shape=19, size=0.5, alpha = 1.0) + 
  geom_smooth(method="lm", se = T, aes(fill=MCQcat_bi), alpha = 0.2) + 
  labs(y="Confidence", x = "logRT (z-score)", color = "MCQ-cat") + 
  theme_minimal() + theme(axis.text=element_text(size=18),axis.title=element_text(size=25))

## 3b) MCQ-feelings

#conduct a hierarchical regression to see if interaction w/ MCQ can explain how much confidence is influenced by standardized log RT
confModel_noMCQ = lmer(conf ~ acc + logRT + acc * logRT + (1 + acc + logRT|subj), data=bigData_clean
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

confModel_wMCQ = lmer(conf ~ MCQ_emo*(acc + logRT + acc * logRT) + (1 + acc + logRT|subj), data=bigData_clean
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wMCQ)
print(summary(confModel_wMCQ))
print(Anova(confModel_wMCQ, type = 3))
coef(summary(confModel_wMCQ)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wMCQ))) 

## check if including MCQ as interaction improved the fit of the model, line 257
anova(confModel_noMCQ,confModel_wMCQ) 

# get median MCQ
MCQ_emo_med <- median(bigData_clean$MCQ_emo)
bigData_clean$MCQemo_bi = cut(as.numeric(bigData_clean$MCQ_emo),breaks=c(0,MCQ_emo_med,max(bigData_clean$MCQ_emo)),include.lowest=T, labels=c("low", "high"))
bigData_clean$MCQemo_bi <- factor(bigData_clean$MCQemo_bi)

#make a nice figure of conf~RTxMCQ interaction 
ggplot(bigData_clean, aes(x=logRT, y=conf, colour=MCQemo_bi)) + 
  geom_count() + 
  scale_color_manual(values=c("salmon", "turquoise3")) +
  geom_point(shape=19, size=0.5, alpha = 1.0) + 
  geom_smooth(method="lm", se = T, aes(fill=MCQemo_bi), alpha = 0.2) + 
  labs(y="Confidence", x = "logRT (z-score)", color = "MCQ-feelings") + 
  theme_minimal() + theme(axis.text=element_text(size=18),axis.title=element_text(size=25))

## 4) EXAMINE INTERACTIONS WITH RAADS

##conduct a hierarchical regression to see if interaction w/ RAADS can explain how much confidence is influenced by standardized log RT
confModel_noRAADS = lmer(conf ~ acc + logRT + acc * logRT + (1 + acc + logRT|subj), data=bigData_clean
                         , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

confModel_wRAADS = lmer(conf ~ RAADS14*(acc + logRT + acc * logRT) + (1 + acc + logRT|subj), data=bigData_clean
                        , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wRAADS)
print(summary(confModel_wRAADS))
print(Anova(confModel_wRAADS, type = 3))
coef(summary(confModel_wRAADS)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wRAADS))) 

## check if including RAADS as interaction improved the fit of the model 
anova(confModel_noRAADS,confModel_wRAADS) 

## 4) EXTRACT CORRECT/ERROR BETAS FOR PLOTTING (NOT USED)

# get median-split on RAADS
bigData_clean$RAADS14 = bigData_clean$RAADS14 - 14 # correct RAADS-14 scores into units of 0-3 on each question to match previous literature 
bigData_clean$RAADS14_menta = bigData_clean$RAADS14_menta - 7 # correct RAADS-14 scores into units of 0-3 on each question to match previous literature 
RAADS_med <- median(bigData_clean$RAADS14)
bigData_clean$RAADSbi = cut(as.numeric(bigData_clean$RAADS14),breaks=c(0,RAADS_med,max(bigData_clean$RAADS14)),include.lowest=T, labels=c("low", "high"))
bigData_clean$RAADSbi <- factor(bigData_clean$RAADSbi)

## distinguish between error/correct & high vs. low autism trials
bigData_err <- bigData_clean[bigData_clean$acc == -0.5, ]
bigData_corr <- bigData_clean[bigData_clean$acc == 0.5, ]
bigData_lASDerr <- bigData_err[bigData_err$RAADSbi == "low", ]
bigData_hASDerr <- bigData_err[bigData_err$RAADSbi == "high", ]
bigData_lASDcorr <- bigData_corr[bigData_corr$RAADSbi == "low", ]
bigData_hASDcorr <- bigData_corr[bigData_corr$RAADSbi == "high", ]

##get the beta coefficients for correct/err
corr_lASD = lmer(conf ~ logRT + (1 + logRT|subj), data=bigData_lASDcorr,
                 control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_lASD)
fix.se <- sqrt(diag(vcov(corr_lASD)))
betas <- c(fix, fix.se)
setwd(dataDir)
write.csv(betas, file = paste('corr_lASD.csv'))

err_lASD = lmer(conf ~ logRT + (1 + logRT|subj), data=bigData_lASDerr,
                 control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_lASD)
fix.se <- sqrt(diag(vcov(err_lASD)))
betas <- c(fix, fix.se)
setwd(dataDir)
write.csv(betas, file = paste('err_lASD.csv'))

corr_hASD = lmer(conf ~ logRT + (1 + logRT|subj), data=bigData_hASDcorr,
                 control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_hASD)
fix.se <- sqrt(diag(vcov(corr_hASD)))
betas <- c(fix, fix.se)
setwd(dataDir)
write.csv(betas, file = paste('corr_hASD.csv'))

err_hASD = lmer(conf ~ logRT + (1 + logRT|subj), data=bigData_hASDerr,
                control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_hASD)
fix.se <- sqrt(diag(vcov(err_hASD)))
betas <- c(fix, fix.se)
setwd(dataDir)
write.csv(betas, file = paste('err_hASD.csv'))


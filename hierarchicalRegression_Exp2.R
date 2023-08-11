## cleans up and analyses autism + comparison data
## elisavanderplasATgmail.com

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
require(ggplot2)
require(plyr)
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe

MCQ_cat = NULL
MCQ_emo = NULL
asdData = NULL 
ctlData = NULL
bigData = NULL
dat = NULL
asdDir = "~/Dropbox/InPrep/Autism/Code/MetaMenta-project-github/Data/Exp2/" 
ctlDir = "~/Dropbox/InPrep/Autism/Code/MetaMenta-project-github/Data/Exp1/"

current_recruitment = "data_exp_27169-v8"
files = c('_task-pf6t', '_task-yzt9')
for (i in 1:length(files)){
  data = read.csv2(paste(asdDir, current_recruitment, '/', current_recruitment, files[i], '.csv', sep = ""), header = T, sep=",",)
  dat = rbind(dat,data)
}

ASD = read.csv2(paste(asdDir, current_recruitment, '/','total_IDs_v8.csv', sep = ""), header = T, sep = ",",)
CTL = read.csv2(paste(asdDir, current_recruitment, '/','selected_comparisons.csv', sep = ""), header = T, sep=",",)
asdIDs = ASD$prolific_ID
ctlIDS = CTL$prolific_ID

for (s in 1:length(asdIDs)){
  subj_dat = dat[dat$Participant.Private.ID==asdIDs[s],] ##load variables for specific subject
  conftask_dat=subj_dat[subj_dat$Task_type=="simpleperceptual",]##select trials from the confidence sub-task
  
  vistrials=conftask_dat[conftask_dat$label=="responsePerceptual",]##select initial binary decision (left/right) trials from the confidence sub-task
  conftrials=conftask_dat[conftask_dat$label=="confidencerating",]##select subsequent confidence rating from the confidence sub-task
  
  logRT = scale(log(vistrials$Reactiontime))
  conf = round(as.numeric(conftrials$confidence_rating)*100)/100
  keypress=vistrials$key_press
  acc=vistrials$correct 
  acc[is.na(acc)]=0##accuracy==1: correct, accuracy==0: wrong
  acc = acc-0.5##accuracy==0.5: correct, accuracy==-0.5: wrong
  
  #have to ecompute objectively correct answer because js script doesn't give that yet
  dir=rep(1, length(acc))
  for (t in 1:length(acc)){
    if (acc[t] ==1 & keypress[t] == 87){
      dir[t]= -1}
    else if (acc[t] == 0 && keypress[t]==69){
      dir[t] = -1}
  }##correct and chose left, dir == -1 (left) or wrong and chose right, dir == -1 (left)
  
  #get all vars behind each other per subject
  subj = rep(s, length(acc))
  group = rep(-0.5, length(acc))
  subData1 = data.frame("subj"=subj, "group"=group, "dir"=dir,"acc"=acc, "conf"=conf, "logRT"=logRT)
  
  #add to larger file 
  asdData = rbind(asdData, subData1)
}

##now do the same for comparisons
dat = NULL
CTLData = NULL
sessions = c('v38', 'v39', 'v40', 'v41', 'v42', 'v43')
#load in all data
for (i in 1:length(sessions)){ 
  currentData = paste('data_exp_12022-',sessions[i],sep="")
  Dir = paste(ctlDir,currentData, '/', sep="")
  files = c('_task-pf6t', '_task-yzt9')
  for (j in 1:length(files)){
    data = read.csv2(paste(Dir, currentData, files[j],'.csv',sep=""), header=T, sep=",", )
    dat = rbind(dat, data)
  }
}
#scale ASD variables
for (s in 1:length(ctlIDS))
{
  subj_dat = dat[dat$Participant.Private.ID==ctlIDS[s],] ##load variables for specific subject
  conftask_dat=subj_dat[subj_dat$Task_type=="simpleperceptual",]##select trials from the confidence sub-task
 
  vistrials=conftask_dat[conftask_dat$label=="responsePerceptual",]##select initial binary decision (left/right) trials from the confidence sub-task
  conftrials=conftask_dat[conftask_dat$label=="confidencerating",]##select subsequent confidence rating from the confidence sub-task

  logRT = scale(log(vistrials$Reactiontime))
  conf = round(as.numeric(conftrials$confidence_rating)*100)/100
  keypress=vistrials$key_press
  acc=vistrials$correct 
  acc[is.na(acc)]=0##accuracy==1: correct, accuracy==0: wrong
  acc = acc-0.5##accuracy==0.5: correct, accuracy==-0.5: wrong
  
  #have to ecompute objectively correct answer because js script doesn't give that yet
  dir=rep(1, length(acc))
  for (t in 1:length(acc)){
    if (acc[t] ==1 & keypress[t] == 87){
      dir[t]= -1}
    else if (acc[t] == 0 && keypress[t]==69){
      dir[t] = -1}
    }##correct and chose left, dir == -1 (left) or wrong and chose right, dir == -1 (left)
    
  #get all vars behind each other per subject
  subj = rep(s, length(acc)) + 40
  group = rep(0.5, length(acc))
  subData2 = data.frame("subj"=subj,"group"=group, "dir"=dir, "acc"=acc, "conf"=conf, "logRT"=logRT)
    
  #add to larger file 
  CTLData = rbind(CTLData, subData2)
}
bigData <- rbind(asdData, CTLData)
# Factors
bigData$subj <- factor(bigData$subj)
bigData$group <- factor(bigData$group, labels=c("ASD", "comparison"))

# skip rows with NaNs
bigData_clean <- na.omit(bigData)

#Conduct a hierarchical regression to see if interaction w/ ASD can explain how much confidence is influenced by standardized log RT
confModel_noASD = lmer(conf ~ acc + logRT + acc * logRT + (1 + acc + logRT|subj), data=bigData_clean
                        , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

confModel_wASD = lmer(conf ~ group*(acc + logRT + acc * logRT) + (1 + acc + logRT|subj), data=bigData_clean
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wASD)
print(summary(confModel_wASD))
print(Anova(confModel_wASD, type = 3))
coef(summary(confModel_wASD)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wASD))) 

## check if including AQ trials as interaction improved the fit of the model 
anova(confModel_noASD,confModel_wASD) 


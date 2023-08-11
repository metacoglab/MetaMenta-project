%% Runs all analyses for the Meta Menta Clinical Study (Experiment2)
%elisavanderplasATgmail.com

clear all; close all; fs = filesep;
findpath = which('DataAnalysis_Exp1.m');
baseDir = fileparts(findpath);
dirData_clinical = [baseDir fs 'Data' fs 'Exp2' fs];
scriptDir = [baseDir 'Analyses' fs];
addpath([baseDir fs 'myfunctions' fs 'HMeta-d' fs 'Matlab' fs]);
addpath([baseDir fs 'myfunctions']);

%load in the prepared data   
load([dirData_clinical 'allClinicalData.mat'])
load([dirData_clinical 'allMetaData.mat'])
load([dirData_clinical 'allDotsData.mat'])

%merge two groups of MetaData together for whole-sample analyses
metaData{3}.nR_S1 = [metaData{1}.nR_S1, metaData{2}.nR_S1];
metaData{3}.nR_S2 = [metaData{1}.nR_S2, metaData{2}.nR_S2];

% Get means and SDs for table (ASD = 1, CTL = 2)
mean_age(1) = nanmean(Data.age(Data.group == 0.5));
sd_age(1) = nanstd(Data.age(Data.group == 0.5));
mean_age(2) = nanmean(Data.age(Data.group ~= 0.5));
sd_age(2) = nanstd(Data.age(Data.group ~= 0.5));
mean_age(1) = nanmean(Data.age(Data.group == 0.5));
sd_age(1) = nanstd(Data.age(Data.group == 0.5));
mean_age(2) = nanmean(Data.age(Data.group ~= 0.5));
sd_age(2) = nanstd(Data.age(Data.group ~= 0.5));
mean_MCQcat(1) = nanmean(Data.MCQ_cat(Data.group == 0.5).*12);
sd_MCQcat(1) = nanstd(Data.MCQ_cat(Data.group == 0.5).*12);
mean_MCQcat(2) = nanmean(Data.MCQ_cat(Data.group ~= 0.5).*12);
sd_MCQcat(2) = nanstd(Data.MCQ_cat(Data.group ~= 0.5).*12);
mean_MCQfeelings(1) = nanmean(Data.MCQ_feelings(Data.group == 0.5).*8);
sd_MCQfeelings(1) = nanstd(Data.MCQ_feelings(Data.group == 0.5).*8);
mean_MCQfeelings(2) = nanmean(Data.MCQ_feelings(Data.group ~= 0.5).*8);
sd_MCQfeelings(2) = nanstd(Data.MCQ_feelings(Data.group ~= 0.5).*8);

%zscore all variables in the table
Data.MCQ_feelings = zscore(Data.MCQ_feelings);
Data.edu = zscore(Data.edu); 
Data.iq = (Data.iq-nanmean(Data.iq))/nanstd(Data.iq); 
Data.age = zscore(Data.age); 

%make a separate dataset without negative Mratios
Data1 = Data(Data.metaR > 0,:);
Data1.metaR = zscore(log(Data1.metaR));

%make a separate dataset without NaN Raads
Data2 = Data(find(~isnan(Data.RAADS)),:);

%independent samples ttest accuracy, line 490
[H,P,CI, STATS] = ttest2(Data.acc(1:40), Data.acc(41:end));
[H,P,KSSTAT] = kstest2(Data.acc(1:40), Data.acc(41:end));

%NB. for the mixed-effect hierarchical regression model, see:
%hierarchicalRegression_Exp2.R

%check menta-ASD impairment effect, line 498
fitlm(Data, 'MCQ_feelings~group+age+gender+edu+IQ')

%% Hypothesis 2
%step 1 - ASD: frequentist linear model w/ covariates
fitlm(Data1, 'metaR~group+age+gender+edu+IQ')

%step 2 - ASD: simultaneous HMeta-d' hierarchical regression
Data4 = Data(find(~isnan(Data.iq)),:);
ids_to_exclude = find(isnan(Data.iq));
j=1;
for i = 1:length(metaData{3}.nR_S1)
    if ~any(i == ids_to_exclude)
        metaData{4}.nR_S1{j} = metaData{3}.nR_S1{i};
        metaData{4}.nR_S2{j} = metaData{3}.nR_S2{i};
        j=j+1;
    end
end
ASD_group = Data4.group == 0.5;

cov_CTL = [Data4.age(~ASD_group)'; Data4.edu(~ASD_group)';Data4.gender(~ASD_group)'; Data4.iq(~ASD_group)'];  
cov_ASD = [Data4.age(ASD_group)'; Data4.edu(ASD_group)';Data4.gender(ASD_group)'; Data4.iq(ASD_group)'];
j=1; k=1;
for i = 1:length(metaData{4}.nR_S1)
    if ASD_group(i) == 1
        metaData_ASD.nR_S1{j} = metaData{4}.nR_S1{i};
        metaData_ASD.nR_S2{j} = metaData{4}.nR_S2{i};
        j=j+1;
    else
        metaData_CTL.nR_S1{k} = metaData{4}.nR_S1{i};
        metaData_CTL.nR_S2{k} = metaData{4}.nR_S2{i};
        k=k+1;
    end
end

FIT2.ASDregr = fit_meta_d_mcmc_regression(metaData_ASD.nR_S1, metaData_ASD.nR_S2, cov_ASD); 
FIT2.CTLregr = fit_meta_d_mcmc_regression(metaData_CTL.nR_S1, metaData_CTL.nR_S2, cov_CTL);
cd(baseDir)
[fig4, fig5] = groupModelfit_checks(FIT2.ASDregr, FIT2.CTLregr,Data.metaR);

% Same group comparison without covariates 
FIT2.ASD = fit_meta_d_mcmc_group(metaData{1}.nR_S1, metaData{1}.nR_S2); 
FIT2.CTL = fit_meta_d_mcmc_group(metaData{2}.nR_S1, metaData{2}.nR_S2);
cd(baseDir)
[fig4, fig5] = groupModelfit_checks(FIT2.ASD, FIT2.CTL,Data.metaR);

%% Runs all analyses for the Meta Menta Study (Experiment 1)
%elisavanderplasATgmail.com

clear all; close all; fs = filesep;
baseDir =  ['~' fs 'Dropbox' fs 'MetaMenta-project' fs];
dirData = [baseDir 'Data' fs 'Exp1' fs];
addpath([baseDir fs 'myfunctions' fs 'HMeta-d' fs 'Matlab' fs]);
addpath([baseDir fs 'myfunctions']);

%load in the prepared data   
load([dirData 'allData.mat'])
load([dirData 'allMetaData.mat'])

%zscore key variables in the table
data.MCQ_feelings = zscore(data.MCQ_feelings);
data.edu = zscore(data.education); 
data.IQ = (data.IQ-nanmean(data.IQ))/nanstd(data.IQ); 
data.age = zscore(data.age); 
data.AQ10 = (data.AQ10-nanmean(data.AQ10))/nanstd(data.AQ10)*-1; %because inv-scored
data.AQ10_C = (data.AQ10_C-nanmean(data.AQ10_C))/nanstd(data.AQ10_C)*-1; %because inv-scored
data.RAADS = (data.RAADS-nanmean(data.RAADS))/nanstd(data.RAADS);
data.RAADS_M = (data.RAADS_MENTA-nanmean(data.RAADS_MENTA))/nanstd(data.RAADS_MENTA);
data.RAADS_NS = (data.RAADS_SOA + data.RAADS_SOR)/2; %take non-social factors of RAADS together
data.RAADS_NS = (data.RAADS_NS-nanmean(data.RAADS_NS))/nanstd(data.RAADS_NS); 

%make a separate dataset without negative Mratios
data1 = data(data.metaR > 0,:);
data1.metaR = zscore(log(data1.metaR));

%make a separate dataset without NaNs in the questionnaires
data2 = data(find(~isnan(data.RAADS)),:); 
metaData(2).nR_S1 = metaData(1).nR_S1;%TO DO remove 241 and 322 (in that order)
metaData(2).nR_S2 = metaData(1).nR_S2;

%% Hypothesis 1
%step 1: simultaneous HMeta-d' hierarchical regression
FIT1 = fit_meta_d_mcmc_regression(metaData(1).nR_S1, metaData(1).nR_S2, data.MCQ_feelings');
cd(baseDir)
[fig1, fig2, fig3] = modelfit_checks(FIT1, data1.metaR, data1.MCQ_feelings, 'Menta');
% compute probability
samples = FIT1.mcmc.samples.mu_beta1 > 0;
p_theta = 1 - ((sum(samples(:) == 0))/(sum(samples(:))));

%step 2: frequentist linear model w/ covariates
fitlm(data1, 'metaR~MCQ_feelings+age+gender+education+IQ')
%NB, for "distinct constructions of confidence in mentalizing" see 'hierarchicalRegression_Exp1.r' script

[fig4] = betaPlot(1); %%plot the coefficients (Figure 3a)

<<<<<<< HEAD
%% replicate menta-AQ effect, line 641-642
=======
%% replicate menta-AQ effect, line 276-277
>>>>>>> d09fc66f474a8bae6dc72fe26517b79d7efa5953
%step 1 - AQ: frequentist linear model w/ covariates
fitlm(data, 'MCQ_cat~AQ10+age+gender+education+IQ')
%step 2 - RAADS: frequentist linear model w/ covariates
fitlm(data, 'MCQ_cat~RAADS+age+gender+education+IQ') %% so focus on the RAADS

%% Hypothesis 2
%step 1 - RAADS: simultaneous HMeta-d' hierarchical regression
FIT2 = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS');
cd(baseDir)
<<<<<<< HEAD
[fig5, fig6, fig7] = modelfit_checks(FIT2,data2.metaR, data2.RAADS, 'RAADS-14');
=======
[fig4, fig5, fig6] = modelfit_checks(FIT2,data2.metaR, data2.RAADS, 'RAADS-14');
>>>>>>> d09fc66f474a8bae6dc72fe26517b79d7efa5953
% compute probability
samples = FIT2.mcmc.samples.mu_beta1 > 0;
p_theta = (sum(samples(:) == 0))/(sum(samples(:)));

%step 2 - RAADS: frequentist linear model w/ covariates
fitlm(data1, 'metaR~RAADS+age+gender+education+IQ')

%NB, for "distinct constructions of confidence in RAADS" see 'hierarchicalRegression_Exp1.r' script

<<<<<<< HEAD
%% Exploratory 3A
%RAADS-Menta: simultaneous HMeta-d' hierarchical regression
FIT4 = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS_M');
cd(baseDir)
[fig8, fig9, fig10] = modelfit_checks(FIT4, data2.metaR, data2.RAADS_M, 'RAADS menta.'); 
% compute probability
samples = FIT4.mcmc.samples.mu_beta1 < 0; 
=======
%% Exploratory 3
%step 1 - AQ10-communication: simultaneous HMeta-d' hierarchical regression
FIT3 = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.AQ10_C');
cd(baseDir)
[fig4, fig5, fig6] = modelfit_checks(FIT3, data2.metaR, data2.AQ10_C, 'AQ10 communication'); 
% compute probability
samples = FIT3.mcmc.samples.mu_beta1 < 0; 
>>>>>>> d09fc66f474a8bae6dc72fe26517b79d7efa5953
p_theta = (sum(samples(:) == 1))/(sum(samples(:))); 

%step 2 - communication: frequentist linear model w/ covariates
fitlm(data1, 'metaR~RAADS_M+age+gender+education+IQ') 

<<<<<<< HEAD
%% Exploratory 3B
%RAADS-nonsocial: simultaneous HMeta-d' hierarchical regression
FIT5 = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS_NS');
cd(baseDir)
[fig4, fig5, fig6] = modelfit_checks(FIT5, data2.metaR, data2.RAADS_NS, 'RAADS non-social.'); 
% compute probability
samples = FIT5.mcmc.samples.mu_beta1 < 0; 
p_theta = (sum(samples(:) == 1))/30000; 
=======
%step 1 - RAADS-Menta: simultaneous HMeta-d' hierarchical regression
FIT4 = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS_M');
cd(baseDir)
[fig4, fig5, fig6] = modelfit_checks(FIT4, data2.metaR, data2.RAADS_M, 'RAADS menta.'); 
% compute probability
samples = FIT4.mcmc.samples.mu_beta1 < 0; 
p_theta = (sum(samples(:) == 1))/(sum(samples(:))); 
>>>>>>> d09fc66f474a8bae6dc72fe26517b79d7efa5953

%step 2 - communication: frequentist linear model w/ covariates
fitlm(data1, 'metaR~RAADS_NS+age+gender+education+IQ') 


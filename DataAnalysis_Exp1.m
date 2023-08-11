%% Runs all analyses for the Meta Menta Study (Experiment 1)
%elisavanderplasATgmail.com

clear all; close all; fs = filesep;
findpath = which('DataAnalysis_Exp1.m');
baseDir = fileparts(findpath);
dirData = [baseDir fs 'Data' fs 'Exp1' fs];
addpath([baseDir fs 'myfunctions' fs 'HMeta-d' fs 'Matlab' fs]);
addpath([baseDir fs 'myfunctions']);

%load in the prepared data   
load([dirData 'allData.mat'])
load([dirData 'allMetaData.mat'])

%% 1) Descriptive statistics on raw data
% correct RAADS-14 scores into units of 0-3 on each question to match previous literature
data.RAADS = data.RAADS - 14;  
data.RAADS_MENTA = data.RAADS_MENTA - 7;
data.RAADS_SOA = data.RAADS_SOA - 4;
data.RAADS_SOR = data.RAADS_SOR - 3;
data.RAADS_NS = (data.RAADS_SOA + data.RAADS_SOR)/2; %take non-social factors of RAADS together

% correct AQ10 scores to negate the double reverse coding error in get_AQ10
data.AQ10 = 10-data.AQ10;

% correct MCQ scores into raw scores for plotting
data.MCQ_cat = data.MCQ_cat.*12;
data.MCQ_feelings = data.MCQ_feelings.*8;

% Get means and SDs for table 
mean_age = nanmean(data.age);
sd_age = nanstd(data.age);
mean_RAADS = nanmean(data.RAADS);
sd_RAADS = nanstd(data.RAADS);
mean_AQ10 = nanmean(data.AQ10);
sd_AQ10 = nanstd(data.AQ10);
mean_MCQcat = nanmean(data.MCQ_cat);
sd_MCQcat = nanstd(data.MCQ_cat);
mean_MCQfeelings = nanmean(data.MCQ_feelings);
sd_MCQfeelings = nanstd(data.MCQ_feelings);
mean_ICAR = nanmean(data.IQ);
sd_ICAR = nanstd(data.IQ);
mean_BCIS = nanmean(data.BCIS);
sd_BCIS = nanstd(data.BCIS);

%zscore key variables in the table for regressions
data_z = data;
data_z.MCQ_feelings = zscore(data.MCQ_feelings);
data_z.MCQ_cat = zscore(data.MCQ_cat);
data_z.edu = zscore(data.education); 
data_z.IQ = (data.IQ-nanmean(data.IQ))/nanstd(data.IQ); 
data_z.age = zscore(data.age); 
data_z.AQ10 = (data.AQ10-nanmean(data.AQ10))/nanstd(data.AQ10);
data_z.RAADS = (data.RAADS-nanmean(data.RAADS))/nanstd(data.RAADS);
data_z.RAADS_M = (data.RAADS_MENTA-nanmean(data.RAADS_MENTA))/nanstd(data.RAADS_MENTA);
data_z.RAADS_NS = (data.RAADS_NS-nanmean(data.RAADS_NS))/nanstd(data.RAADS_NS);

%make separate datasets without outlier negative Mratios and take log
data1 = data(data.metaR > 0.01,:);
data1.metaR = log(data1.metaR);
data1_z = data_z(data_z.metaR > 0.01,:);
data1_z.metaR = log(data1_z.metaR);

%make a separate dataset without NaNs in the questionnaires
data2 = data_z(find(~isnan(data.RAADS)),:); 
metaData(1) = metaData;
% adjust the metacognitive data for input into the regression to remove
% these individuals with no Qs data
ids_to_exclude = find(isnan(data.RAADS));
j=1;
for i = 1:length(metaData(1).nR_S1)
    if ~any(i == ids_to_exclude)
        metaData(2).nR_S1{j} = metaData(1).nR_S1{i};
        metaData(2).nR_S2{j} = metaData(1).nR_S2{i};
        j=j+1;
    end
end

%% 2) correlations between autism scores and mentalising scores 
[rho_AQ_RAADS,pval_AQ_RAADS] = corr(data2.AQ10, data2.RAADS, 'Type', 'Spearman');
[rho_MCQcat_RAADS,pval_MCQcat_RAADS] = corr(data2.MCQ_cat, data2.RAADS, 'Type', 'Spearman');
[rho_MCQfeelings_RAADS,pval_MCQfeelings_RAADS] = corr(data2.MCQ_feelings, data2.RAADS, 'Type', 'Spearman');
[rho_MCQcat_AQ10,pval_MCQcat_AQ10] = corr(data2.MCQ_cat, data2.AQ10, 'Type', 'Spearman');
[rho_MCQfeelings_AQ10 ,pval_MCQfeelings_AQ10] = corr(data2.MCQ_feelings, data2.AQ10, 'Type', 'Spearman');

% linear models controlling for covariates relating MCQ_cat / feelings to both AQ10
% and RAADS
fitlm(data_z, 'MCQ_cat~AQ10+age+gender+education+IQ')
fitlm(data_z, 'MCQ_cat~RAADS+age+gender+education+IQ')
fitlm(data_z, 'MCQ_feelings~AQ10+age+gender+education+IQ')
fitlm(data_z, 'MCQ_feelings~RAADS+age+gender+education+IQ')

%% 3) Relationships with mentalising

%% 3a) Relationship between metacognitive efficiency and MCQ_categorisation 
%step 1: simultaneous HMeta-d' hierarchical regression
fit_MCQcat = fit_meta_d_mcmc_regression(metaData(1).nR_S1, metaData(1).nR_S2, data_z.MCQ_cat');
cd(baseDir)
[fig1, fig2] = modelfit_checks(fit_MCQcat, 'MCQ-cat');
% compute probability
p_theta_mcqcat = sum(fit_MCQcat.mcmc.samples.mu_beta1(:) > 0)./(length(fit_MCQcat.mcmc.samples.mu_beta1(:)));

%step 2: frequentist linear model w/ covariates
fitlm(data1_z, 'metaR~MCQ_cat+age+gender+education+IQ')

% Plot single-subject relationship
figure;
scatter(data1.metaR, data1.MCQ_cat, 70,'Marker', 'o', 'MarkerFaceColor',[0.5, 0.5, 0.5],'LineWidth',2);
hLine = refline; 
hLine.Color = 'k'; 
hLine.LineWidth = 5;
xlabel('log Mratio', 'FontSize', 34);
ylabel(['MCQ-cat'], 'FontSize',34);
set(gca, 'FontSize',24)
set(gcf, 'color', 'w');

%% 3b) Relationship between metacognitive efficiency and MCQ_feelings
%% Hypothesis 1
%step 1: simultaneous HMeta-d' hierarchical regression
fit_MCQfeelings = fit_meta_d_mcmc_regression(metaData(1).nR_S1, metaData(1).nR_S2, data_z.MCQ_feelings');
cd(baseDir)
[fig3, fig4] = modelfit_checks(fit_MCQfeelings, 'MCQ-feelings');
% compute probability
p_theta_mcqfeelings = sum(fit_MCQfeelings.mcmc.samples.mu_beta1(:) > 0)./(length(fit_MCQfeelings.mcmc.samples.mu_beta1(:)));

%step 2: frequentist linear model w/ covariates
fitlm(data1_z, 'metaR~MCQ_feelings+age+gender+education+IQ')
%NB, for "distinct constructions of confidence in mentalizing" see 'hierarchicalRegression_Exp1.r' script

% Plot single-subject relationship
figure;
scatter(data1.metaR, data1.MCQ_feelings, 70,'Marker', 'o', 'MarkerFaceColor',[0.5, 0.5, 0.5],'LineWidth',2);
hLine = refline; 
hLine.Color = 'k'; 
hLine.LineWidth = 5;
xlabel('log Mratio', 'FontSize', 34);
ylabel(['MCQ-feelings'], 'FontSize',34);
set(gca, 'FontSize',24)
set(gcf, 'color', 'w');

%% 4) Relationship between metacognitive efficiency and RAADS

%step 1 - RAADS: simultaneous HMeta-d' hierarchical regression
fit_RAADS = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS');
cd(baseDir)
[fig5, fig6] = modelfit_checks(fit_RAADS, 'RAADS-14');
% compute probability (negative, so use less than)
p_theta_RAADS = sum(fit_RAADS.mcmc.samples.mu_beta1(:) < 0)./(length(fit_RAADS.mcmc.samples.mu_beta1(:)));

%step 2 - RAADS: frequentist linear model w/ covariates
fitlm(data1_z, 'metaR~RAADS+age+gender+education+IQ')

%% 4a) Relationship between metacognitive efficiency and RAADS-M

fit_RAADS_M = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS_MENTA');
cd(baseDir)
[fig7, fig8] = modelfit_checks(fit_RAADS_M, 'RAADS-14 mentalising');
% compute probability (negative, so use less than)
p_theta_RAADS_M = sum(fit_RAADS_M.mcmc.samples.mu_beta1(:) < 0)./(length(fit_RAADS_M.mcmc.samples.mu_beta1(:)));

%step 2 - RAADS: frequentist linear model w/ covariates
fitlm(data1_z, 'metaR~RAADS_MENTA+age+gender+education+IQ')

%% 4b) Relationship between metacognitive efficiency and RAADS-NS

fit_RAADS_NS = fit_meta_d_mcmc_regression(metaData(2).nR_S1, metaData(2).nR_S2, data2.RAADS_NS');
cd(baseDir)
[fig9, fig10] = modelfit_checks(fit_RAADS_NS, 'RAADS-14 non-mentalising');
% compute probability (negative, so use less than)
p_theta_RAADS_NS = sum(fit_RAADS_NS.mcmc.samples.mu_beta1(:) < 0)./(length(fit_RAADS_NS.mcmc.samples.mu_beta1(:)));

%step 2 - RAADS: frequentist linear model w/ covariates
fitlm(data1_z, 'metaR~RAADS_NS+age+gender+education+IQ')
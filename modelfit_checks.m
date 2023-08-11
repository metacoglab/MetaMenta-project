function [fig1, fig2] = modelfit_checks(MODEL, name_var1)

%% takes as input a HMeta-d regression model fit, as well as the group level metacognition (var1) and the dependent variable (var2)
%% returns the trace plot and posterior density histogram plot (fig1), and a correlation plot between var1 and var2 (fig2), 

%% fig 2 is the same as the metad-group-visualise.m in the HMeta-d' toolbox

color_TOM = {[0.129, 0.050, 0.568], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

fig1 = figure;
subplot(2,2,1)
plot(exp(MODEL.mcmc.samples.mu_logMratio')) 
title(['Rhat: ' num2str(MODEL.mcmc.Rhat.mu_logMratio)], 'fontsize', 12);
xlabel('Sample', 'fontsize', 14);
ylabel('\mu_{meta-d''/d''}', 'fontsize', 14);
box off

subplot(2,2,2)
histogram(exp(MODEL.mcmc.samples.mu_logMratio(:)),500, 'facecolor', [0.6, 0.6, 0.6], 'edgecolor', [0.6, 0.6, 0.6], 'facealpha', 0.4);
title('Group-level Mratio', 'fontsize', 12); 
HDI = calc_HDI(exp(MODEL.mcmc.samples.mu_logMratio(:)));
% text(HDI(1), 235, num2str(round(HDI,2)+0.01), 'FontSize', 14); 
xline(mean(exp(MODEL.mcmc.samples.mu_logMratio(:))), '-', 'color', 'b', 'linewidth', 2)
xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', 14);
xlabel('\mu_{meta-d''/d''}', 'FontSize',14);
set(gca, 'XLim', [0.65 0.85], 'YLim', [0 300]);
box off

subplot(2,2,3)
hold all
plot(MODEL.mcmc.samples.mu_beta1') 
title(['Rhat: ' num2str(MODEL.mcmc.Rhat.mu_beta1)], 'fontsize', 12);
xlabel('Sample', 'fontsize', 14);
ylabel(name_var1, 'fontsize', 14);
box off

subplot(2,2,4)
histogram(MODEL.mcmc.samples.mu_beta1(:),500,'facecolor', [0.6, 0.6, 0.6], 'edgecolor', [0.6, 0.6, 0.6], 'facealpha', 0.4);
title(['Beta of regression with ' name_var1], 'fontsize', 12); 
HDI = calc_HDI(MODEL.mcmc.samples.mu_beta1(:));
xline(0, '-', 'color', 'k', 'linewidth', 2)
ylabel('No. of samples', 'FontSize', 14);
xlabel(name_var1, 'FontSize', 14);
xlim([-0.05 0.1])
xline(0, '-', 'color', 'k', 'linewidth', 1)
set(gca, 'XLim', [-0.05 0.15], 'YLim', [0 300]);
set(gcf, 'color', 'w')
box off

fig2 = figure; 
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.2 0.2 0.5 0.5]);
set(gcf, 'color', 'w')

Nsub = length(MODEL.d1);
ts = tinv([0.05/2,  1-0.05/2],Nsub-1);

if any(isnan(MODEL.obs_FAR2_rS1(:))) || any(isnan(MODEL.obs_HR2_rS1(:))) || any(isnan(MODEL.obs_FAR2_rS2(:))) || any(isnan(MODEL.obs_HR2_rS2(:)))
    warning('One or more subjects have NaN entries for observed confidence rating counts; these will be omitted from the plot')
end

mean_obs_FAR2_rS1 = nanmean(MODEL.obs_FAR2_rS1);
mean_obs_HR2_rS1 = nanmean(MODEL.obs_HR2_rS1);
mean_obs_FAR2_rS2 = nanmean(MODEL.obs_FAR2_rS2);
mean_obs_HR2_rS2 = nanmean(MODEL.obs_HR2_rS2);

CI_obs_FAR2_rS1(1,:) = ts(1).*(nanstd(MODEL.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS1(2,:) = ts(2).*(nanstd(MODEL.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(1,:) = ts(1).*(nanstd(MODEL.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(2,:) = ts(2).*(nanstd(MODEL.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS2(1,:) = ts(1).*(nanstd(MODEL.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_FAR2_rS2(2,:) = ts(2).*(nanstd(MODEL.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(1,:) = ts(1).*(nanstd(MODEL.obs_HR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(2,:) = ts(2).*(nanstd(MODEL.obs_HR2_rS2)./sqrt(Nsub));

mean_est_FAR2_rS1 = nanmean(MODEL.est_FAR2_rS1);
mean_est_HR2_rS1 = nanmean(MODEL.est_HR2_rS1);
mean_est_FAR2_rS2 = nanmean(MODEL.est_FAR2_rS2);
mean_est_HR2_rS2 = nanmean(MODEL.est_HR2_rS2);

CI_est_FAR2_rS1(1,:) = ts(1).*(nanstd(MODEL.est_FAR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS1(2,:) = ts(2).*(nanstd(MODEL.est_FAR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(1,:) = ts(1).*(nanstd(MODEL.est_HR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(2,:) = ts(2).*(nanstd(MODEL.est_HR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS2(1,:) = ts(1).*(nanstd(MODEL.est_FAR2_rS2)./sqrt(Nsub));
CI_est_FAR2_rS2(2,:) = ts(2).*(nanstd(MODEL.est_FAR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(1,:) = ts(1).*(nanstd(MODEL.est_HR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(2,:) = ts(2).*(nanstd(MODEL.est_HR2_rS2)./sqrt(Nsub));

%% Observed and expected type 2 ROCs for S1 and S2 responses
subplot(1,2,1);
errorbar([1 mean_obs_FAR2_rS1 0], [1 mean_obs_HR2_rS1 0], [0 CI_obs_HR2_rS1(1,:) 0], [0 CI_obs_HR2_rS1(2,:) 0], ...
    [0 CI_obs_FAR2_rS1(1,:) 0], [0 CI_obs_FAR2_rS1(2,:) 0], 'ko-','linewidth',1.5,'markersize', 12);
hold on
errorbar([1 mean_est_FAR2_rS1 0], [1 mean_est_HR2_rS1 0], [0 CI_est_HR2_rS1(1,:) 0], [0 CI_est_HR2_rS1(2,:) 0], ...
    [0 CI_est_FAR2_rS1(1,:) 0], [0 CI_est_FAR2_rS1(2,:) 0], 'd-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2');
xlabel('FAR2');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off
legend('Data', 'Model', 'Location', 'SouthEast')
legend boxoff
title('Response = S1')

subplot(1,2,2);
errorbar([1 mean_obs_FAR2_rS2 0], [1 mean_obs_HR2_rS2 0], [0 CI_obs_HR2_rS2(1,:) 0], [0 CI_obs_HR2_rS2(2,:) 0], ...
    [0 CI_obs_FAR2_rS2(1,:) 0], [0 CI_obs_FAR2_rS2(2,:) 0], 'ko-','linewidth',1.5,'markersize', 12);
hold on
errorbar([1 mean_est_FAR2_rS2 0], [1 mean_est_HR2_rS2 0], [0 CI_est_HR2_rS2(1,:) 0], [0 CI_est_HR2_rS2(2,:) 0], ...
    [0 CI_est_FAR2_rS2(1,:) 0], [0 CI_est_FAR2_rS2(2,:) 0], 'd-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2');
xlabel('FAR2');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off
legend('Data', 'Model', 'Location', 'SouthEast')
legend boxoff
title('Response = S2')


end

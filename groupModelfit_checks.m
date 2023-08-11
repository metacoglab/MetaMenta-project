function [fig1, fig2] = groupModelfit_checks(MODEL1, MODEL2, var1)
%% takes as input two independent group model fit, as well as the group level metacognition (var1) and the dependent variable (var2)
%% returns the trace plot and posterior density histogram plot (fig1), and a correlation plot between var1 and var2 (fig2)
color_TOM = {[0.427, 0.298, 0.803], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

sampleDiff = MODEL1.mcmc.samples.mu_logMratio - MODEL2.mcmc.samples.mu_logMratio;
hdi = calc_HDI(sampleDiff(:));
hdi_asd = calc_HDI(exp(MODEL1.mcmc.samples.mu_logMratio))
hdi_cnt = calc_HDI(exp(MODEL2.mcmc.samples.mu_logMratio))

fprintf(('\n Mratio group values = %.2f and %.2f'), exp(MODEL1.mu_logMratio), exp(MODEL2.mu_logMratio));
fprintf(['\n Estimated difference in Mratio between groups: ', num2str(exp(MODEL1.mu_logMratio) - exp(MODEL2.mu_logMratio))])
fprintf(['\n HDI on difference in log(Mratio): ', num2str(hdi) '\n\n'])

% compute probability of difference
temp = sampleDiff < 0;
p_theta = (sum(temp(:) == 1))/30000;
fprintf(('\n Probability that sample difference is lower than zero: %.2f'), p_theta); 

fig1 = figure;
hold on 
for i = 1:40
    r = -0.065 + (0.065+0.065).*rand(1,1);
    scat1 = plot([1+r], var1(i),  'o', 'MarkerSize', 9, 'MarkerFaceColor', color_TOM{1}, 'MarkerEdgeColor', color_TOM{1});
end
for i = 1:40
    r = -0.065 + (0.065+0.065).*rand(1,1);
    scat2 = plot([2+r], var1(i+40),  'o', 'MarkerSize', 9, 'MarkerFaceColor',color_TOM{2}, 'MarkerEdgeColor', color_TOM{2});
end
b1 = bar([1 2],[nanmean(var1(1:40)), nanmean(var1(41:end))], 'FaceColor', 'w', 'BarWidth', 0.8, 'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.4);
errorbar([1 2],[nanmean(var1(1:40)), nanmean(var1(41:end))],[nanmean(var1(1:40))./sqrt(40), nanstd(var1(41:end))./sqrt(40)], '.', 'Color', [0 0 0], 'LineWidth', 4)
set(gca, 'YLim', [0 1], 'XLim', [0  3], 'XTick', [1,2], 'XTickLabel', {'ASD', 'CNTRL'}, 'FontSize', 20);
ylabel([{'Standardized individual fit'};{'Metacognitive efficiency'}], 'FontSize', 20);
xlabel('Clinical group', 'FontSize', 20);
set(gcf, 'color', 'w')

[lgd, handle] = legend([scat1, scat2], {'autism', 'comparison'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',15);

fig2 = figure;
hold all
title('95% HDI', 'fontsize', 20); 
HDI = calc_HDI(exp(MODEL1.mcmc.samples.mu_logMratio(:)));
leg1=xline(HDI(1), '--', 'color', color_TOM{1}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{1}, 'linewidth', 2)
histogram(exp(MODEL1.mcmc.samples.mu_logMratio(:)),500, 'facecolor', color_TOM{1}, 'facealpha', 0.4, 'edgecolor', color_TOM{1}, 'edgealpha', 0.4);
histogram(exp(MODEL2.mcmc.samples.mu_logMratio(:)),500,'facecolor', color_TOM{2}, 'facealpha', 0.4, 'edgecolor', color_TOM{2}, 'edgealpha', 0.4);
HDI = calc_HDI(exp(MODEL2.mcmc.samples.mu_logMratio(:)));
leg2 = xline(HDI(1), '--', 'color', color_TOM{2}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{2}, 'linewidth', 2)
xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', 18);
set(gca, 'XLim', [-0.05 0.15], 'YLim', [0 300], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',18);
xlim([-0.2 1.2])
set(gcf, 'color', 'w')

[lgd, handle] = legend([leg1, leg2], {'autism', 'comparison'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',15);
hold off

Nsub = length(MODEL2.d1);
ts = tinv([0.05/2,  1-0.05/2],Nsub-1);

if any(isnan(MODEL2.obs_FAR2_rS1(:))) || any(isnan(MODEL2.obs_HR2_rS1(:))) || any(isnan(MODEL2.obs_FAR2_rS2(:))) || any(isnan(MODEL2.obs_HR2_rS2(:)))
    warning('One or more subjects have NaN entries for observed confidence rating counts; these will be omitted from the plot')
end

mean_obs_FAR2_rS1 = nanmean(MODEL2.obs_FAR2_rS1);
mean_obs_HR2_rS1 = nanmean(MODEL2.obs_HR2_rS1);
mean_obs_FAR2_rS2 = nanmean(MODEL2.obs_FAR2_rS2);
mean_obs_HR2_rS2 = nanmean(MODEL2.obs_HR2_rS2);

CI_obs_FAR2_rS1(1,:) = ts(1).*(nanstd(MODEL2.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS1(2,:) = ts(2).*(nanstd(MODEL2.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(1,:) = ts(1).*(nanstd(MODEL2.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(2,:) = ts(2).*(nanstd(MODEL2.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS2(1,:) = ts(1).*(nanstd(MODEL2.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_FAR2_rS2(2,:) = ts(2).*(nanstd(MODEL2.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(1,:) = ts(1).*(nanstd(MODEL2.obs_HR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(2,:) = ts(2).*(nanstd(MODEL2.obs_HR2_rS2)./sqrt(Nsub));

mean_est_FAR2_rS1 = nanmean(MODEL2.est_FAR2_rS1);
mean_est_HR2_rS1 = nanmean(MODEL2.est_HR2_rS1);
mean_est_FAR2_rS2 = nanmean(MODEL2.est_FAR2_rS2);
mean_est_HR2_rS2 = nanmean(MODEL2.est_HR2_rS2);

CI_est_FAR2_rS1(1,:) = ts(1).*(nanstd(MODEL2.est_FAR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS1(2,:) = ts(2).*(nanstd(MODEL2.est_FAR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(1,:) = ts(1).*(nanstd(MODEL2.est_HR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(2,:) = ts(2).*(nanstd(MODEL2.est_HR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS2(1,:) = ts(1).*(nanstd(MODEL2.est_FAR2_rS2)./sqrt(Nsub));
CI_est_FAR2_rS2(2,:) = ts(2).*(nanstd(MODEL2.est_FAR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(1,:) = ts(1).*(nanstd(MODEL2.est_HR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(2,:) = ts(2).*(nanstd(MODEL2.est_HR2_rS2)./sqrt(Nsub));

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
legend('Data', 'MODEL2', 'Location', 'SouthEast')
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
legend('Data', 'MODEL2', 'Location', 'SouthEast')
legend boxoff
title('Response = S2')


end

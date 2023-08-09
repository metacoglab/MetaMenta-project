function [h01] = betaPlot(exp)

%% Plot the beta coefficients from R (see: HierarchicalRegression_Exp1.R) to figure "h01"
%% receives as input whether needs to plot the beta plot of experiment 1 or 2, as exp = 1/2
% EVDP 2019 elisa.plas.18@ucl.ac.uk

fs = filesep;
h01 = figure;
c.corr =  [0.082, 0.615, 0.835];
c.err = [0.835, 0.250, 0.082];
color_TOM = {[0.427, 0.298, 0.803], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

if exp == 1
    baseDir =  ['~' fs 'Dropbox' fs 'MetaMenta-project' fs 'Data' fs 'Exp1' fs];
    
    for acc = 1:2
        for ASD = 1:2
            
            correct = {'corr_'; 'err_'};
            ASDtrait = {'hASD'; 'lASD'};
            suffix = [correct{acc} ASDtrait{ASD}];
            
            %% Load betas
            cd(baseDir)
            datafile = [suffix '.csv'];
            dat{ASD, acc} = readtable(datafile);
        end
    end
    
    ms = 22;
    axis_text = 26;
    axis_nr = 18;
    hold all
    plot(1,  [dat{1,1}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.corr, 'Linewidth', 3, 'MarkerEdgeColor', c.corr);
    plot(1,  [dat{1,2}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
    x1 = plot(2,  [dat{2,1}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.corr,'Linewidth', 3, 'MarkerEdgeColor', c.corr);
    x2= plot(2,  [dat{2,2}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
    errorbar(1:2, [dat{1,1}{2,2} dat{2,1}{2,2}], [dat{1,1}{4,2} dat{2,1}{4,2}], '.', 'Color', c.corr*0.2, 'LineWidth', 2);
    errorbar(1:2, [dat{1,2}{2,2} dat{2,2}{2,2}], [dat{1,2}{4,2} dat{2,2}{4,2}], '.', 'Color', c.err*0.2 , 'LineWidth', 2);
    
    hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline
    set(gca, 'XLim', [0 3], 'XTick', 1:2, 'YLim', [-0.55 0.05],'YTick', -0.55:0.1:0.5, 'FontSize',axis_nr, 'XTickLabels', {'High RAADS', 'Low RAADS'}, 'Fontsize', 18);
    ylabel([{'logRT impact'};{'on confidence (a.u.)'}], 'FontSize', axis_text);
    set(gcf, 'color', 'w')
    
    
[lgd, handle] = legend([x1, x2],{'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',1)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);
    
else
    
    baseDir =  ['~' fs 'Dropbox' fs 'MetaMenta-project' fs 'Data' fs 'Exp2' fs];
    filename = 'regression_betas_'; %'HierarchicalRegression_exp2.r'
    
    for acc = 1:2
        for ASD = 1:2
            correct = {'_corr'; '_err'};
            clinical = {'ASD'; 'CTL'};
            suffix = [clinical{ASD} correct{acc}];
            
            %% Load betas
            cd(baseDir)
            datafile = [filename suffix '.csv'];
            dat{ASD, acc} = readtable(datafile);
        end
    end
    
    ms = 23;
    axis_text = 26;
    axis_nr = 18;
    hold all
    plot(1,  [dat{1,1}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', color_TOM{1},'Linewidth', 3, 'MarkerEdgeColor', 'k');
    plot(1,  [dat{1,2}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', color_TOM{1}, 'Linewidth', 3, 'MarkerEdgeColor','k');
    plot(2,  [dat{2,1}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', color_TOM{2},'Linewidth', 3, 'MarkerEdgeColor', 'k');
    plot(2,  [dat{2,2}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', color_TOM{2}, 'Linewidth', 3, 'MarkerEdgeColor','k');
    errorbar(1:2, [dat{1,1}{2,2} dat{2,1}{2,2}], [dat{1,1}{2,2}./sqrt(80) dat{2,1}{2,2}./sqrt(80)], '.', 'Color','k', 'LineWidth', 2);
    errorbar(1:2, [dat{1,2}{2,2} dat{2,2}{2,2}], [dat{1,2}{2,2}./sqrt(80) dat{2,2}{2,2}./sqrt(80)], '.', 'Color', 'k' , 'LineWidth', 2);
    
    fakeline1 = plot(1,  [dat{1,1}{2,2}]+3000 , 'd', 'MarkerSize', 16, 'MarkerFaceColor', 'k','Linewidth', 1, 'MarkerEdgeColor', 'k');
    fakeline2 = plot(1,  [dat{1,2}{2,2}]+30000 , 'o', 'MarkerSize', 16, 'MarkerFaceColor', 'k', 'Linewidth', 1, 'MarkerEdgeColor','k');
    
    hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline
     set(gca, 'XLim', [0 3], 'XTick', 1:2, 'YLim', [-0.55 0.05],'YTick', -0.55:0.1:0.5, 'FontSize',axis_nr, 'XTickLabels', {'ASD', 'CNTRL'}, 'Fontsize', 18);
   ylabel([{'logRT impact'};{'on confidence (a.u.)'}], 'FontSize', axis_text);
    set(gcf, 'color', 'w')
    
    [lgd, handle] = legend([fakeline1, fakeline2],{'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',1)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

end
end

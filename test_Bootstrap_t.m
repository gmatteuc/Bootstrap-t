% \\\\\\\\\\\\\\\\\\\\\\ test_Bootstrap_t \\\\\\\\\\\\\\\\\\\\\\

% TEST BOOTSTRAP-t Test the accompanying code providing a flexible parallel
% implementation of the Bootstrap-t procedure to compute confidence intervals 
% (CI) following the procedure described by the book:
% "Tibshirani, R. J., & Efron, B. (1993).
% An introduction to the bootstrap.
% Monographs on statistics and applied probability, 57(1)"
%
% AUTHOR: Giulio Matteucci
% DATE: 21/04/2023

%% create test dataset to test Bootstrap-t functions

% generate normal data
Npoints=5000;
Xinput1_mu = 0.33;
Xinput1_sig = 0.33;
Xinput1 = normrnd(Xinput1_mu,Xinput1_sig,Npoints,1);
Yinput1 = normrnd(0,1,Npoints,1);
Xinput2_mu = -0.33;
Xinput2_sig = 0.33;
Xinput2 = normrnd(Xinput2_mu,Xinput2_sig,Npoints,1);
Yinput2 = normrnd(0,1,Npoints,1);
quadrant = 1;
takeabs = 0;

% set bootstrap function input data
InputData = {Xinput1, Yinput1, Xinput2, Yinput2, quadrant, takeabs};

% get mean differences
mean_diff_output = get_mean_difference_Bootstrap_t(InputData);
% get ttest2 confidence interval
[ttest2_h,ttest2_p,ttest2_ci,ttest2_stats] = ttest2(InputData{1},InputData{3});

% plot the histograms of input data
f1=figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2, 1, 1);
histogram(InputData{1},50,'EdgeAlpha',0,'FaceColor',[1,0,0].*0.75);
hold on;
histogram(InputData{3},50,'EdgeAlpha',0,'FaceColor',[0,1,0].*0.75);
xlabel('value');
ylabel('frequency');
title(['input data distributions ( Xinput1mu = ',num2str(Xinput1_mu),' - Xinput2mu = ',num2str(Xinput2_mu),' )']);
legend('Xinput1', 'Xinput2');
% plot the mean difference as a vertical line
y_limits = ylim;
y_limits = [y_limits(1),1.1*y_limits(2)];
line([mean_diff_output{2}, mean_diff_output{2}], y_limits, 'Color', [1,0,0].*0.75, 'LineWidth', 2, 'LineStyle', '--');
line([mean_diff_output{3}, mean_diff_output{3}], y_limits, 'Color', [0,1,0].*0.75, 'LineWidth', 2, 'LineStyle', '--');
hold off;
% display the mean difference in the command window
disp(['Mean Difference: ', num2str(mean_diff_output{1})]);
xlim([-2.5,2.5]);
ylim(y_limits);
grid on;
set(gca,'fontsize',12)

%% perform Bootstrap-t results

% initialize output structures
orig_estimate=cell(1,3); %#ok<*PREALL>
orig_estimate_lCI=cell(1,3);
orig_estimate_uCI=cell(1,3);
orig_estimate_lSE=cell(1,3);
orig_estimate_uSE=cell(1,3);
new_estimate=cell(1,3);
new_estimate_lCI=cell(1,3);
new_estimate_uCI=cell(1,3);
new_estimate_lSE=cell(1,3);
new_estimate_uSE=cell(1,3);
toc1=NaN; %#ok<*NASGU>
toc2=NaN;

% set bootstrap function input data
Binp = {Xinput1, Yinput1, Xinput2, Yinput2, quadrant, takeabs};
% set bootstrap function input parameters
confidence = 0.95;
B = 1500;
N = 25;
Bfunc = @get_mean_difference_Bootstrap_t;
Bdim = [1, 1, 2, 2, 0, 0];
Brdim = [1, 1, 1, 1, 0, 0];

tic1=tic;
% call the (original) for-based function
[orig_estimate, orig_estimate_lCI, orig_estimate_uCI, orig_estimate_lSE, orig_estimate_uSE] =...
    get_Bootstrap_t_ci_serial(Bfunc, Binp, Bdim, Brdim, confidence, B, N, 1, 1);
% maesure elapsed time (eta)
toc1=toc(tic1);

tic2=tic;
% call the (new) parfor-based function
[new_estimate, new_estimate_lCI, new_estimate_uCI, new_estimate_lSE, new_estimate_uSE] =...
    get_Bootstrap_t_ci_parallel(Bfunc, Binp, Bdim, Brdim, confidence, B, N, 1, 1);
% maesure elapsed time (eta)
toc2=toc(tic2);

%% inspect Bootstrap-t results

% set output names
outputnames={'mean difference','mean 1','mean 2'};

% display result comparison
disp(['----------------------------------']);
disp(['orig eta: ', num2str(toc1), ' s new eta: ', num2str(toc2),' s']);
for i = 1:length(orig_estimate)
    disp(['output ', num2str(i), ' (', outputnames{i} ') --------------']);
    if i==1
        disp(['orig estimate: ', num2str(orig_estimate{i}), ' new estimate: ', num2str(new_estimate{i}),' ground truth: ',num2str(Xinput1_mu-Xinput2_mu)]);
        disp(['orig lCI: ', num2str(orig_estimate_lCI{i}), ' new lCI: ', num2str(new_estimate_lCI{i}), ' ttest2 lCI: ', num2str(ttest2_ci(1))]);
        disp(['orig uCI: ', num2str(orig_estimate_uCI{i}), ' new uCI: ', num2str(new_estimate_uCI{i}), ' ttest2 uCI: ', num2str(ttest2_ci(2))]);
    else
        disp(['orig estimate: ', num2str(orig_estimate{i}), ' new estimate: ', num2str(new_estimate{i})]);
        disp(['orig lCI: ', num2str(orig_estimate_lCI{i}), ' new lCI: ', num2str(new_estimate_lCI{i})]);
        disp(['orig uCI: ', num2str(orig_estimate_uCI{i}), ' new uCI: ', num2str(new_estimate_uCI{i})]);
    end
end
disp(['----------------------------------']);

% prepare for plotting result comparison
estimates = [orig_estimate{1}, new_estimate{1}, mean_diff_output{1}];
lower_CI = [orig_estimate_lCI{1}, new_estimate_lCI{1}, ttest2_ci(1)];
upper_CI = [orig_estimate_uCI{1}, new_estimate_uCI{1}, ttest2_ci(2)];
lower_error = estimates - lower_CI;
upper_error = upper_CI - estimates;
errors = [lower_error; upper_error];
x = 1:length(estimates);
xlimtouse=[0,4];
% plot result comparison
f2=figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
scatter(x,estimates,75,'Markerfacecolor',[0,0,0],'Markeredgecolor',[0,0,0]);
errorbar(x,estimates,lower_error,upper_error,'LineStyle','none','Linewidth',2,'Color',[0,0,0]);
pl=plot([xlimtouse(1),xlimtouse(2)],[Xinput1_mu-Xinput2_mu,Xinput1_mu-Xinput2_mu],'linewidth',1,'color',[0.5,0.5,0.5]);
hold off;
xticks(x);
xtickangle(45);
xticklabels({'Bootstrap-t (serial)', 'Bootstrap-t (parallel)', 't-test'});
ylabel('estimate value');
title(['Bootstrap-t and parametric estimates with 95% CI - ',outputnames{1}]);
ylim(gca,[0.4,1]);
xlim(gca,xlimtouse)
grid on;
legend(pl,'ground truth');
set(gca,'fontsize',12)

% --------------------------------------------------------------------
% Validation (performed on MATLAB R2019b):
% Example results of a test run comparing the CI from an unpaired ttest with the
% bootstrap-estimated CI (Npoints=5000, B=10000, N=100, sig=0.33, mu=+/-0.33)
% with "get_boostrap_parallel_fast_t_ci_BN2p" function:
% --------------------------------------------------------------------
% new eta: 109.231 s
% mean difference --------------
% Bootstrap-t estimate: 0.66774 - true value = 0.66000
% Bootstrap-t lCI: 0.65496 - ttest2 lCI: 0.65495
% Bootstrap-t uCI: 0.68065 - ttest2 uCI: 0.68053
% --------------------------------------------------------------------
% Disclaimer:
% This software is provided "as-is" without any warranty, express or implied.
% The author(s) and/or contributors shall not be held liable for any damages,
% including but not limited to, direct, indirect, incidental, consequential,
% or other losses, resulting from the use of, or inability to use, this software.
% By using this software, you acknowledge that you understand and agree to this disclaimer.
% Use at your own risk.
% --------------------------------------------------------------------
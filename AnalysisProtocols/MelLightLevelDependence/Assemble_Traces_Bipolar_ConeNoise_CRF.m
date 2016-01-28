clear all; close all; clc;

contrastLevels = [9 18 36];

%%%%% LMS+ %%%%%
%% Load CSV files
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sLMSDirectedPositivePulseConeNoiseCRF_';
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_09pct.csv']));
    M_LMSPos_09(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_18pct.csv']));
    M_LMSPos_18(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+_36pct.csv']));
    M_LMSPos_36(:, s) = tmp(:, 2);
end

avgLMSPos_09 = mean(M_LMSPos_09, 2); semLMSPos_09 = std(M_LMSPos_09, [], 2)/sqrt(size(M_LMSPos_09, 2));
avgLMSPos_18 = mean(M_LMSPos_18, 2); semLMSPos_18 = std(M_LMSPos_18, [], 2)/sqrt(size(M_LMSPos_18, 2));
avgLMSPos_36 = mean(M_LMSPos_36, 2); semLMSPos_36 = std(M_LMSPos_36, [], 2)/sqrt(size(M_LMSPos_36, 2));

subplot(1, 3, 1);
shadedErrorBar(M_t(1:600), 100*avgLMSPos_09, 100*semLMSPos_09); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+ 9%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 2);
shadedErrorBar(M_t(1:600), 100*avgLMSPos_18, 100*semLMSPos_18); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+ 18%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 3);
shadedErrorBar(M_t(1:600), 100*avgLMSPos_36, 100*semLMSPos_36); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+ 36%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

set(gcf, 'PaperPosition', [0 0 7 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [7 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'LMSPos_CRF_TimeSeries.pdf'), 'pdf');
close(gcf)

%% Summary statistics
% Min pupil size
subplot(1, 3, 1);
minLMSPos_09 = min(M_LMSPos_09(100:200, :));
minLMSPos_18 = min(M_LMSPos_18(100:200, :));
minLMSPos_36 = min(M_LMSPos_36(100:200, :));
avgMinLMSPos_09 = mean(minLMSPos_09);
avgMinLMSPos_18 = mean(minLMSPos_18);
avgMinLMSPos_36 = mean(minLMSPos_36);
semMinLMSPos_09 = std(minLMSPos_09, [], 2)/sqrt(size(minLMSPos_09, 2));
semMinLMSPos_18 = std(minLMSPos_18, [], 2)/sqrt(size(minLMSPos_36, 2));
semMinLMSPos_36 = std(minLMSPos_36, [], 2)/sqrt(size(minLMSPos_18, 2));

plot(log2(contrastLevels), 100*[avgMinLMSPos_09 avgMinLMSPos_18 avgMinLMSPos_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinLMSPos_09 avgMinLMSPos_18 avgMinLMSPos_36], 100*[semMinLMSPos_09 semMinLMSPos_18 semMinLMSPos_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-0.35 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('Peak pupillary constriction [%]');

% AUC
subplot(1, 3, 2);
aucLMSPos_09 = trapz(M_t(100:200), M_LMSPos_09(100:200, :));
aucLMSPos_18 = trapz(M_t(100:200), M_LMSPos_18(100:200, :));
aucLMSPos_36 = trapz(M_t(100:200), M_LMSPos_36(100:200, :));
avgMinLMSPos_09 = mean(aucLMSPos_09);
avgMinLMSPos_18 = mean(aucLMSPos_18);
avgMinLMSPos_36 = mean(aucLMSPos_36);
semMinLMSPos_09 = std(aucLMSPos_09, [], 2)/sqrt(size(aucLMSPos_09, 2));
semMinLMSPos_18 = std(aucLMSPos_18, [], 2)/sqrt(size(aucLMSPos_36, 2));
semMinLMSPos_36 = std(aucLMSPos_36, [], 2)/sqrt(size(aucLMSPos_18, 2));

plot(log2(contrastLevels), 100*[avgMinLMSPos_09 avgMinLMSPos_18 avgMinLMSPos_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinLMSPos_09 avgMinLMSPos_18 avgMinLMSPos_36], 100*[semMinLMSPos_09 semMinLMSPos_18 semMinLMSPos_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-1 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('AUC');
close(gcf)


%%%%% LMS- %%%%%
%% Load CSV files
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sLMSDirectedNegativePulseConeNoiseCRF_';
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS-_09pct.csv']));
    M_LMSNeg_09(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS-_18pct.csv']));
    M_LMSNeg_18(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS-_36pct.csv']));
    M_LMSNeg_36(:, s) = tmp(:, 2);
end

avgLMSNeg_09 = mean(M_LMSNeg_09, 2); semLMSNeg_09 = std(M_LMSNeg_09, [], 2)/sqrt(size(M_LMSNeg_09, 2));
avgLMSNeg_18 = mean(M_LMSNeg_18, 2); semLMSNeg_18 = std(M_LMSNeg_18, [], 2)/sqrt(size(M_LMSNeg_18, 2));
avgLMSNeg_36 = mean(M_LMSNeg_36, 2); semLMSNeg_36 = std(M_LMSNeg_36, [], 2)/sqrt(size(M_LMSNeg_36, 2));

subplot(1, 3, 1);
shadedErrorBar(M_t(1:600), 100*avgLMSNeg_09, 100*semLMSNeg_09); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS- 9%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 2);
shadedErrorBar(M_t(1:600), 100*avgLMSNeg_18, 100*semLMSNeg_18); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS- 18%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 3);
shadedErrorBar(M_t(1:600), 100*avgLMSNeg_36, 100*semLMSNeg_36); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS- 36%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

set(gcf, 'PaperPosition', [0 0 7 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [7 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'LMSNeg_CRF_TimeSeries.pdf'), 'pdf');
close(gcf)

%% Summary statistics
% Min pupil size
subplot(1, 3, 1);
minLMSNeg_09 = min(M_LMSNeg_09(100:200, :));
minLMSNeg_18 = min(M_LMSNeg_18(100:200, :));
minLMSNeg_36 = min(M_LMSNeg_36(100:200, :));
avgMinLMSNeg_09 = mean(minLMSNeg_09);
avgMinLMSNeg_18 = mean(minLMSNeg_18);
avgMinLMSNeg_36 = mean(minLMSNeg_36);
semMinLMSNeg_09 = std(minLMSNeg_09, [], 2)/sqrt(size(minLMSNeg_09, 2));
semMinLMSNeg_18 = std(minLMSNeg_18, [], 2)/sqrt(size(minLMSNeg_36, 2));
semMinLMSNeg_36 = std(minLMSNeg_36, [], 2)/sqrt(size(minLMSNeg_18, 2));

plot(log2(contrastLevels), 100*[avgMinLMSNeg_09 avgMinLMSNeg_18 avgMinLMSNeg_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinLMSNeg_09 avgMinLMSNeg_18 avgMinLMSNeg_36], 100*[semMinLMSNeg_09 semMinLMSNeg_18 semMinLMSNeg_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-0.35 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('Peak pupillary constriction [%]');

% AUC
subplot(1, 3, 2);
aucLMSNeg_09 = trapz(M_t(100:200), M_LMSNeg_09(100:200, :));
aucLMSNeg_18 = trapz(M_t(100:200), M_LMSNeg_18(100:200, :));
aucLMSNeg_36 = trapz(M_t(100:200), M_LMSNeg_36(100:200, :));
avgMinLMSNeg_09 = mean(aucLMSNeg_09);
avgMinLMSNeg_18 = mean(aucLMSNeg_18);
avgMinLMSNeg_36 = mean(aucLMSNeg_36);
semMinLMSNeg_09 = std(aucLMSNeg_09, [], 2)/sqrt(size(aucLMSNeg_09, 2));
semMinLMSNeg_18 = std(aucLMSNeg_18, [], 2)/sqrt(size(aucLMSNeg_36, 2));
semMinLMSNeg_36 = std(aucLMSNeg_36, [], 2)/sqrt(size(aucLMSNeg_18, 2));

plot(log2(contrastLevels), 100*[avgMinLMSNeg_09 avgMinLMSNeg_18 avgMinLMSNeg_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinLMSNeg_09 avgMinLMSNeg_18 avgMinLMSNeg_36], 100*[semMinLMSNeg_09 semMinLMSNeg_18 semMinLMSNeg_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-1 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('AUC');
close(gcf)


%%%%% Mel- %%%%%
%% Load CSV files
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnoreNegativePulseConeNoiseCRF_';
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel-_09pct.csv']));
    M_MelNeg_09(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel-_18pct.csv']));
    M_MelNeg_18(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel-_36pct.csv']));
    M_MelNeg_36(:, s) = tmp(:, 2);
end

avgMelNeg_09 = mean(M_MelNeg_09, 2); semMelNeg_09 = std(M_MelNeg_09, [], 2)/sqrt(size(M_MelNeg_09, 2));
avgMelNeg_18 = mean(M_MelNeg_18, 2); semMelNeg_18 = std(M_MelNeg_18, [], 2)/sqrt(size(M_MelNeg_18, 2));
avgMelNeg_36 = mean(M_MelNeg_36, 2); semMelNeg_36 = std(M_MelNeg_36, [], 2)/sqrt(size(M_MelNeg_36, 2));

subplot(1, 3, 1);
shadedErrorBar(M_t(1:600), 100*avgMelNeg_09, 100*semMelNeg_09); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel- 9%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 2);
shadedErrorBar(M_t(1:600), 100*avgMelNeg_18, 100*semMelNeg_18); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel- 18%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 3);
shadedErrorBar(M_t(1:600), 100*avgMelNeg_36, 100*semMelNeg_36); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel- 36%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

set(gcf, 'PaperPosition', [0 0 7 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [7 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'MelNeg_CRF_TimeSeries.pdf'), 'pdf');
close(gcf)

%% Summary statistics
% Min pupil size
subplot(1, 3, 1);
minMelNeg_09 = min(M_MelNeg_09(100:200, :));
minMelNeg_18 = min(M_MelNeg_18(100:200, :));
minMelNeg_36 = min(M_MelNeg_36(100:200, :));
avgMinMelNeg_09 = mean(minMelNeg_09);
avgMinMelNeg_18 = mean(minMelNeg_18);
avgMinMelNeg_36 = mean(minMelNeg_36);
semMinMelNeg_09 = std(minMelNeg_09, [], 2)/sqrt(size(minMelNeg_09, 2));
semMinMelNeg_18 = std(minMelNeg_18, [], 2)/sqrt(size(minMelNeg_36, 2));
semMinMelNeg_36 = std(minMelNeg_36, [], 2)/sqrt(size(minMelNeg_18, 2));

plot(log2(contrastLevels), 100*[avgMinMelNeg_09 avgMinMelNeg_18 avgMinMelNeg_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinMelNeg_09 avgMinMelNeg_18 avgMinMelNeg_36], 100*[semMinMelNeg_09 semMinMelNeg_18 semMinMelNeg_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-0.35 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('Peak pupillary constriction [%]');

% AUC
subplot(1, 3, 2);
aucMelNeg_09 = trapz(M_t(100:200), M_MelNeg_09(100:200, :));
aucMelNeg_18 = trapz(M_t(100:200), M_MelNeg_18(100:200, :));
aucMelNeg_36 = trapz(M_t(100:200), M_MelNeg_36(100:200, :));
avgMinMelNeg_09 = mean(aucMelNeg_09);
avgMinMelNeg_18 = mean(aucMelNeg_18);
avgMinMelNeg_36 = mean(aucMelNeg_36);
semMinMelNeg_09 = std(aucMelNeg_09, [], 2)/sqrt(size(aucMelNeg_09, 2));
semMinMelNeg_18 = std(aucMelNeg_18, [], 2)/sqrt(size(aucMelNeg_36, 2));
semMinMelNeg_36 = std(aucMelNeg_36, [], 2)/sqrt(size(aucMelNeg_18, 2));

plot(log2(contrastLevels), 100*[avgMinMelNeg_09 avgMinMelNeg_18 avgMinMelNeg_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinMelNeg_09 avgMinMelNeg_18 avgMinMelNeg_36], 100*[semMinMelNeg_09 semMinMelNeg_18 semMinMelNeg_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-1 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('AUC');
close(gcf)




%%%%% Mel+ %%%%%
%% Load CSV files
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnorePositivePulseConeNoiseCRF_';
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_09pct.csv']));
    M_MelPos_09(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_18pct.csv']));
    M_MelPos_18(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+_36pct.csv']));
    M_MelPos_36(:, s) = tmp(:, 2);
end

avgMelPos_09 = mean(M_MelPos_09, 2); semMelPos_09 = std(M_MelPos_09, [], 2)/sqrt(size(M_MelPos_09, 2));
avgMelPos_18 = mean(M_MelPos_18, 2); semMelPos_18 = std(M_MelPos_18, [], 2)/sqrt(size(M_MelPos_18, 2));
avgMelPos_36 = mean(M_MelPos_36, 2); semMelPos_36 = std(M_MelPos_36, [], 2)/sqrt(size(M_MelPos_36, 2));

subplot(1, 3, 1);
shadedErrorBar(M_t(1:600), 100*avgMelPos_09, 100*semMelPos_09); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+ 9%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 2);
shadedErrorBar(M_t(1:600), 100*avgMelPos_18, 100*semMelPos_18); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+ 18%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(1, 3, 3);
shadedErrorBar(M_t(1:600), 100*avgMelPos_36, 100*semMelPos_36); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+ 36%', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

set(gcf, 'PaperPosition', [0 0 7 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [7 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'MelPos_CRF_TimeSeries.pdf'), 'pdf');
close(gcf)

%% Summary statistics
% Min pupil size
subplot(1, 3, 1);
minMelPos_09 = min(M_MelPos_09(100:200, :));
minMelPos_18 = min(M_MelPos_18(100:200, :));
minMelPos_36 = min(M_MelPos_36(100:200, :));
avgMinMelPos_09 = mean(minMelPos_09);
avgMinMelPos_18 = mean(minMelPos_18);
avgMinMelPos_36 = mean(minMelPos_36);
semMinMelPos_09 = std(minMelPos_09, [], 2)/sqrt(size(minMelPos_09, 2));
semMinMelPos_18 = std(minMelPos_18, [], 2)/sqrt(size(minMelPos_36, 2));
semMinMelPos_36 = std(minMelPos_36, [], 2)/sqrt(size(minMelPos_18, 2));

plot(log2(contrastLevels), 100*[avgMinMelPos_09 avgMinMelPos_18 avgMinMelPos_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinMelPos_09 avgMinMelPos_18 avgMinMelPos_36], 100*[semMinMelPos_09 semMinMelPos_18 semMinMelPos_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-0.35 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('Peak pupillary constriction [%]');

% AUC
subplot(1, 3, 2);
aucMelPos_09 = trapz(M_t(100:200), M_MelPos_09(100:200, :));
aucMelPos_18 = trapz(M_t(100:200), M_MelPos_18(100:200, :));
aucMelPos_36 = trapz(M_t(100:200), M_MelPos_36(100:200, :));
avgMinMelPos_09 = mean(aucMelPos_09);
avgMinMelPos_18 = mean(aucMelPos_18);
avgMinMelPos_36 = mean(aucMelPos_36);
semMinMelPos_09 = std(aucMelPos_09, [], 2)/sqrt(size(aucMelPos_09, 2));
semMinMelPos_18 = std(aucMelPos_18, [], 2)/sqrt(size(aucMelPos_36, 2));
semMinMelPos_36 = std(aucMelPos_36, [], 2)/sqrt(size(aucMelPos_18, 2));

plot(log2(contrastLevels), 100*[avgMinMelPos_09 avgMinMelPos_18 avgMinMelPos_36], '-ok', 'MarkerFaceColor', 'k'); hold on;
errorbar(log2(contrastLevels), 100*[avgMinMelPos_09 avgMinMelPos_18 avgMinMelPos_36], 100*[semMinMelPos_09 semMinMelPos_18 semMinMelPos_36], 'Color', 'k', 'LineStyle', 'none');
xlim([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)]);
ylim(100*[-1 0.05]);
pbaspect([1 1 1]);
plot([log2(min(contrastLevels)/2) log2(max(contrastLevels)*2)], [0 0], '--', 'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', log2(contrastLevels));
set(gca, 'XTickLabel', {'9%', '18%', '36%'});
xlabel('Contrast'); ylabel('AUC');
close(gcf)
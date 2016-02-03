clear all; close all; clc;
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sBipolarPulse_';


%% Load CSV files
ndVal = 'ND15';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMSPos(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);

        tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS-.csv']));
    M_LMSNeg(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_MelPos(:, s) = tmp(:, 2);
    
            tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel-.csv']));
    M_MelNeg(:, s) = tmp(:, 2);
end

avgMelPos = mean(M_MelPos, 2); semMelPos = std(M_MelPos, [], 2)/sqrt(size(M_MelPos, 2));
avgLMSPos = mean(M_LMSPos, 2); semLMSPos = std(M_LMSPos, [], 2)/sqrt(size(M_LMSPos, 2));
avgMelNeg = mean(M_MelNeg, 2); semMelNeg = std(M_MelNeg, [], 2)/sqrt(size(M_MelNeg, 2));
avgLMSNeg = mean(M_LMSNeg, 2); semLMSNeg = std(M_LMSNeg, [], 2)/sqrt(size(M_LMSNeg, 2));

subplot(4, 2, 1);
shadedErrorBar(M_t(1:600), 100*avgLMSPos, 100*semLMSPos); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2, 3);
shadedErrorBar(M_t(1:600), 100*avgLMSNeg, 100*semLMSNeg); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS-', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2, 5);
shadedErrorBar(M_t(1:600), 100*avgMelPos, 100*semMelPos); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+', ndVal});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2, 7);
shadedErrorBar(M_t(1:600), 100*avgMelNeg, 100*semMelNeg); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel-', ndVal});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');
ylabel('Pupil amplitude [\Delta%]');

%% Load CSV files
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMSPos(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);

        tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS-.csv']));
    M_LMSNeg(:, s) = tmp(:, 2);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_MelPos(:, s) = tmp(:, 2);
    
            tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel-.csv']));
    M_MelNeg(:, s) = tmp(:, 2);
end

avgMelPos = mean(M_MelPos, 2); semMelPos = std(M_MelPos, [], 2)/sqrt(size(M_MelPos, 2));
avgLMSPos = mean(M_LMSPos, 2); semLMSPos = std(M_LMSPos, [], 2)/sqrt(size(M_LMSPos, 2));
avgMelNeg = mean(M_MelNeg, 2); semMelNeg = std(M_MelNeg, [], 2)/sqrt(size(M_MelNeg, 2));
avgLMSNeg = mean(M_LMSNeg, 2); semLMSNeg = std(M_LMSNeg, [], 2)/sqrt(size(M_LMSNeg, 2));

subplot(4, 2, 2);
shadedErrorBar(M_t(1:600), 100*avgLMSPos, 100*semLMSPos); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
%ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2, 4);
shadedErrorBar(M_t(1:600), 100*avgLMSNeg, 100*semLMSNeg); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS-', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
%ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2,6);
shadedErrorBar(M_t(1:600), 100*avgMelPos, 100*semMelPos); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+', ndVal});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
%ylabel('Pupil amplitude [\Delta%]');

subplot(4, 2, 8);
shadedErrorBar(M_t(1:600), 100*avgMelNeg, 100*semMelNeg); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel-', ndVal});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');
%ylabel('Pupil amplitude [\Delta%]');


set(gcf, 'PaperPosition', [0 0 3.5 7]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [3.5 7]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'Bipolar_ND10_ND15.pdf'), 'pdf');
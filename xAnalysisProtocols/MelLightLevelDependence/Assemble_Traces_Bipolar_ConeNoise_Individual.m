clear all; close all; clc;
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sBipolarPulseConeNoise_';

theCols = [72 123 107 ; 238 159 96 ; 120 52 54 ]/255;

%% Load CSV files
ndVal = 'ND10';
theSubjects = {'G100815Ax' 'J100715Rx' 'M100915Sx'};
%theSubjects = {'J100715Rx' 'M100915Sx'};
nSubjects = length(theSubjects);

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

subplot(2, 2, 1);
plot(M_t(1:600), 100*avgLMSPos, '-k'); hold on;
for s = 1:nSubjects
    plot(M_t(1:600), 100*M_LMSPos(:, s), '-', 'Color', theCols(s, :), 'LineWidth', 1.0); hold on;
end
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');

subplot(2, 2, 2);
plot(M_t(1:600), 100*avgLMSNeg, '-k'); hold on;
for s = 1:nSubjects
    plot(M_t(1:600), 100*M_LMSNeg(:, s), '-', 'Color', theCols(s, :), 'LineWidth', 1.0); hold on;
end
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS-', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
%ylabel('Pupil amplitude [\Delta%]');

subplot(2, 2,3);
plot(M_t(1:600), 100*avgMelPos, '-k'); hold on;
for s = 1:nSubjects
    plot(M_t(1:600), 100*M_MelPos(:, s), '-', 'Color', theCols(s, :), 'LineWidth', 1.0); hold on;
end
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');
ylabel('Pupil amplitude [\Delta%]');

subplot(2, 2, 4);
plot(M_t(1:600), 100*avgMelNeg, '-k'); hold on;
for s = 1:nSubjects
    plot(M_t(1:600), 100*M_MelNeg(:, s), '-', 'Color', theCols(s, :), 'LineWidth', 1.0); hold on;
end
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel-'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');
%ylabel('Pupil amplitude [\Delta%]');


set(gcf, 'PaperPosition', [0 0 3.5 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [3.5 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'Bipolar_ConeNoise_ND10_Individual.pdf'), 'pdf');
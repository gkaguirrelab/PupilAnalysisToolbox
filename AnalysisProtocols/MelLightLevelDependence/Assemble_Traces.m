clear all; close all; clc;
basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sPulse_';


%% Load CSV files
ndVal = 'ND20';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_Mel(:, s) = tmp(:, 2);
end

avgMel = mean(M_Mel, 2); semMel = std(M_Mel, [], 2)/sqrt(size(M_Mel, 2));
avgLMS = mean(M_LMS, 2); semLMS = std(M_LMS, [], 2)/sqrt(size(M_LMS, 2));

subplot(2, 5, 1);
shadedErrorBar(M_t(1:600), 100*avgLMS, 100*semLMS); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
ylabel('Pupil amplitude [\Delta%]');


subplot(2, 5, 6);
shadedErrorBar(M_t(1:600), 100*avgMel, 100*semMel); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');
ylabel('Pupil amplitude [\Delta%]');

Min_ND20 = min(M_LMS);

%% Load CSV files
ndVal = 'ND15';
theSubjects = {'G100515Ax' 'J100615Rx' 'M100615Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_Mel(:, s) = tmp(:, 2);
end

avgMel = mean(M_Mel, 2); semMel = std(M_Mel, [], 2)/sqrt(size(M_Mel, 2));
avgLMS = mean(M_LMS, 2); semLMS = std(M_LMS, [], 2)/sqrt(size(M_LMS, 2));

subplot(2, 5, 2);
shadedErrorBar(M_t(1:600), 100*avgLMS, 100*semLMS); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);


subplot(2, 5, 7);
shadedErrorBar(M_t(1:600), 100*avgMel, 100*semMel); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');

Min_ND15 = min(M_LMS);

%% Load CSV files
ndVal = 'ND10';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_Mel(:, s) = tmp(:, 2);
end

avgMel = mean(M_Mel, 2); semMel = std(M_Mel, [], 2)/sqrt(size(M_Mel, 2));
avgLMS = mean(M_LMS, 2); semLMS = std(M_LMS, [], 2)/sqrt(size(M_LMS, 2));

subplot(2, 5, 3);
shadedErrorBar(M_t(1:600), 100*avgLMS, 100*semLMS); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);


subplot(2, 5, 8);
shadedErrorBar(M_t(1:600), 100*avgMel, 100*semMel); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');

Min_ND10 = min(M_LMS);

%% Load CSV files
ndVal = 'ND05';
theSubjects = {'G100515Ax' 'J100615Rx' 'M100615Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_Mel(:, s) = tmp(:, 2);
end

avgMel = mean(M_Mel, 2); semMel = std(M_Mel, [], 2)/sqrt(size(M_Mel, 2));
avgLMS = mean(M_LMS, 2); semLMS = std(M_LMS, [], 2)/sqrt(size(M_LMS, 2));

subplot(2, 5, 4);
shadedErrorBar(M_t(1:600), 100*avgLMS, 100*semLMS); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);


subplot(2, 5, 9);
shadedErrorBar(M_t(1:600), 100*avgMel, 100*semMel); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');

Min_ND05 = min(M_LMS);

%% Load CSV files
ndVal = 'ND00';
theSubjects = {'G092815Ax' 'J092915Rx' 'M092515Sx'};

for s = 1:length(theSubjects)
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_LMS+.csv']));
    M_LMS(:, s) = tmp(:, 2);
    M_t = tmp(:, 1);
    
    tmp = csvread(fullfile([basePath ndVal], [theSubjects{s} ndVal '_PupilPulseData_Mel+.csv']));
    M_Mel(:, s) = tmp(:, 2);
end

avgMel = mean(M_Mel, 2); semMel = std(M_Mel, [], 2)/sqrt(size(M_Mel, 2));
avgLMS = mean(M_LMS, 2); semLMS = std(M_LMS, [], 2)/sqrt(size(M_LMS, 2));

subplot(2, 5, 5);
shadedErrorBar(M_t(1:600), 100*avgLMS, 100*semLMS); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'LMS+', ndVal});
plot([5 10], 100*[0.20 0.20], 'r', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);


subplot(2, 5, 10);
shadedErrorBar(M_t(1:600), 100*avgMel, 100*semMel); hold on;
plot([M_t(1) M_t(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title({'Mel+'});
plot([5 10], 100*[0.20 0.20], 'b', 'LineWidth', 3);
ylim(100*[-0.35 0.35]);
xlim([0 25]);
xlabel('Time [s]');

Min_ND00 = min(M_LMS);

set(gcf, 'PaperPosition', [0 0 7 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [7 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'Unipolar_AllND.pdf'), 'pdf');


%% Plot the minimum value for L+M+S
min_LMS = [Min_ND20 ; Min_ND15 ; Min_ND10 ; Min_ND05; Min_ND00];
mean_min_LMS = mean(min_LMS, 2);
sem_min_LMS = std(min_LMS, [], 2);
errorbar([1:5], -100*mean_min_LMS, 100*sem_min_LMS, '-k'); hold on;
plot([1:5], -100*mean_min_LMS, 'ok', 'MarkerFaceColor', 'k');
ylim(100*[0 0.35]);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');
title('L+M+S Minimum');
set(gca, 'XTick', [1:5]); set(gca, 'XTickLabel', {'ND2.0' 'ND1.5' 'ND1.0' 'ND0.5' 'ND0.0'});
xlabel('Light level'); ylabel('Max. pupil constriction [%]');
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [3.5 3.5]); %Set the paper to have width 15 and height 6.
saveDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/MelLightLevelDependence/Plots';
saveas(gcf, fullfile(saveDir, 'LMSMinimum_AllND.png'), 'png');
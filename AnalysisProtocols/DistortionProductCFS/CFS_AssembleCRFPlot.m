function crfPlot
%% crfPlot

%Create line graph measuring log contrats vs. amplitude
clear all; close all;
%% Set Subjects
subjects = {'s001' 's002'};
outputDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/DistortionProductCFS/Plots';
inputDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/analysis/results/PupillometryLightFluxDistortion4HzCRF/';
hold all;
theMarkerColors = [0.667 0.224 0.224 ; 0.41 0.424 0.424];

for i = 1:length(subjects)
    %% Open CSV and load subject data
    dataPath = strcat(inputDir, subjects(i), '-results.csv');
    data = csvread(char(dataPath),1,2);
    contrasts = [2 4 8 16 32 64];
    amplitudes = data(1:length(contrasts),3)*100;
    amplitudeErrors = data(1:length(contrasts),5)*100;
    
    % Fit exponential
    f1 = fit(contrasts',amplitudes,'exp1')
    xfit = contrasts(1):contrasts(end);
    y = feval(f1, xfit);
    
    %% Plot Results
    errorbar(log2(contrasts),amplitudes,amplitudeErrors, 'LineStyle', 'none', 'Color', theMarkerColors(i, :)); hold on;
    plot(log2(contrasts), amplitudes, 'o', 'Color', theMarkerColors(i, :), 'MarkerFaceColor', theMarkerColors(i, :));
    plot(log2(xfit), y, 'Color', theMarkerColors(i, :));
    
    set(gca,'XTick', log2([contrasts]));
    set(gca,'XTickLabel', [contrasts]);
    xlabel('Contrast [%]');
    ylabel('Amplitude [\Delta%]');
    ylim([0 2.6]);
    xlim(log2([1 128]));
    pbaspect([1 1 1]); box off;
    
    set(gcf, 'PaperPosition', [0 0 4 4])
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
    saveas(gcf, fullfile(outputDir, ['CFS_AssembleCRFPlot_' subjects{i} '.pdf']), 'pdf');
    close(gcf);
end

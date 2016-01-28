function ttfPlot
%% ttfPlot

%Create line graph measuring log carrier frequency vs. change in amplitude
clear all; close all;
%% Set Subjects
subjects = {'s001' 's002' 's003' 's004'};
carrierFreqs = [4 8 16 32 64 128];
totalAmplitudes = zeros(length(carrierFreqs),1);
totalAmplitudeErrors = zeros(length(carrierFreqs),1);
outputDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/DistortionProductCFS/Plots';
inputDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/analysis/results/PupillometryDistortionProductTTF/';

theMarkerColors = [0.667 0.224 0.224 ; 0.667 0.424 0.424 ; 0.133 0.4 0.4 ; 0.176 0.533 0.176];

hold all;
for i = 1:length(subjects)
    %% Open CSV and load subject data
    dataPath = strcat(inputDir, subjects(i), '-results.csv');
    data = csvread(char(dataPath),1,2);
    amplitudes = data(1:length(carrierFreqs),3)*100;
    amplitudeErrors = data(1:length(carrierFreqs),5)*100;
    
    totalAmplitudes = totalAmplitudes + amplitudes;
    totalAmplitudeErrors = totalAmplitudeErrors + amplitudeErrors;
    
    %% Plot Results
    shadedErrorBar(log2(carrierFreqs),amplitudes,amplitudeErrors); hold on;
    plot(log2(carrierFreqs),amplitudes, '-o', 'Color', theMarkerColors(i, :), 'MarkerFaceColor', theMarkerColors(i, :));
    
    
    set(gca,'XTick', log2([carrierFreqs]));
    set(gca,'XTickLabel', [carrierFreqs]);
    
    title({subjects{i} '\pm1SEM (bootstrapped across runs)'});
    xlabel('Carrier frequency [Hz]');
    ylabel('Amplitude [\Delta%]');
    xlim(log2([2 256]));
    ylim([0 2.6]);
    pbaspect([1 1 1]); box off;
    
    set(gcf, 'PaperPosition', [0 0 4 4])
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
    saveas(gcf, fullfile(outputDir, ['CFS_AssembleTTFPlot_' subjects{i} '.pdf']), 'pdf');
    close(gcf);
    
end

%% Plot Group Average
avgAmplitudes = totalAmplitudes./length(subjects);
avgAmplitudeErrors = totalAmplitudeErrors./sqrt(length(subjects));

shadedErrorBar(log2(carrierFreqs),avgAmplitudes,avgAmplitudeErrors); hold on;
plot(log2(carrierFreqs), avgAmplitudes, 'ok', 'MarkerFaceColor', 'k');
set(gca,'XTick', log2([carrierFreqs]));
set(gca,'XTickLabel', [carrierFreqs]);
title({'Group average' '\pm1SEM (across subjects)'});

xlabel('Carrier frequency [Hz]');
ylabel('Amplitude [\Delta%]');
xlim(log2([2 256]));
ylim([0 2.6]);
pbaspect([1 1 1]); box off;

set(gcf, 'PaperPosition', [0 0 4 4])
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(outputDir, 'CFS_AssembleTTFPlot_GrpAverage.pdf'), 'pdf');
close(gcf);
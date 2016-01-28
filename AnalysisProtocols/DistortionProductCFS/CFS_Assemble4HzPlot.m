function distortionProductIntro
%% distortionProductIntro

% Bar graph showing results of subjects' amplitude response at 4 Hz to
% introduce concept of distortion product

subjects = {'s001' 's002' 's003' 's004'};
hold all;
amplitudes = zeros(1,length(subjects)+1);
amplitudeErrors = zeros(1,length(subjects)+1);
outputDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/DistortionProductCFS/Plots';
inputDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/analysis/results/PupillometryDistortionProductTTF/';

for i = 1:length(subjects)
    %% Open CSV and load subject data
    dataPath = strcat(inputDir, subjects(i), '-results.csv');
    data = csvread(char(dataPath),1,2);
    amplitudes(i) = data(1,3)*100;
    amplitudeErrors(i) = data(1,5)*100;
end
amplitudes(length(subjects)+1) = mean(amplitudes(1:length(subjects))); %Group Average
amplitudeErrors(length(subjects)+1) = sum(amplitudeErrors(1:length(subjects)))/sqrt(length(subjects));
graph = bar(amplitudes);
hold on;
    ylim([0 2.6]);
set(gca, 'XTick', 1:(length(subjects)+1));
set(gca, 'XTickLabel', [subjects(1:4) 'Group']);
xlabel('Subject');
ylabel('Amplitude [\Delta%]');
title('Amplitude at 4 Hz AM Flicker');
errorbar(amplitudes, amplitudeErrors, '.k');

box off;
pbaspect([1 1 1]);
set(gcf, 'PaperPosition', [0 0 4 4])
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(outputDir, ['CFS_Assemble4HzPlot.pdf']), 'pdf');
close(gcf);
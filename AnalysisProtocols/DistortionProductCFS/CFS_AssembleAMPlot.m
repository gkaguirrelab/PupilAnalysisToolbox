function cfsBarGraph
%% cfsBarGraph

%Create bar graph showing effects of CFS
clear all; close all;
%% Set Subjects
subjects = {'s001' 's002' 's005'};
modulations = {'Background' '0.5/16 Hz' '0.5/4 Hz' '0.5 Hz'};
totalAmplitudes = zeros(length(modulations)*2,1);
totalAmplitudeErrors = zeros(length(modulations)*2,1); %used for computing group average
outputDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/DistortionProductCFS/Plots';
inputDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/analysis/results/PupillometryCFSAttenuationBattery/';

for i = 1:length(subjects)
    
    %% Open CSV and load subject data
    dataPath = strcat(inputDir, subjects(i),'-results.csv');
    data = csvread(char(dataPath),1,2);
    amplitudes = data(1:(length(modulations)*2),3)*100;
    amplitudeErrors = data(1:(length(modulations)*2),5)*100;
    totalAmplitudes = totalAmplitudes + amplitudes;
    totalAmplitudeErrors = totalAmplitudeErrors + amplitudeErrors;
    
    staticBackground = amplitudes(1);
    cfsBackground = amplitudes(2);
    %% Plot Results
    for j = 2:length(modulations)
        cfsData(j-1,:) = [amplitudes((2*j-1):(2*j))]';
        cfsErrorData(j-1,:) = [amplitudeErrors((2*j-1):(2*j))]';
    end
    graph = bar(cfsData);
    hold on;
    set(gca,'XTick', 1:length(modulations)-1);
    set(gca,'XTickLabel', modulations(2:4));
    xlabel('Modulation');
    set(gca, 'YLim', [0 6]);
    ylabel('Amplitude [\Delta%]');
    
    title(subjects(i));
    plot(xlim, [staticBackground staticBackground], '--k');
    plot(xlim, [cfsBackground cfsBackground], '--g');
    hold on;
    set(graph,'BarWidth',1);    % The bars will now touch each other
    numgroups = size(cfsData, 1);
    numbars = size(cfsData, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for j = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, cfsData(:,j), cfsErrorData(:,j), 'k', 'linestyle', 'none');
    end
    legend('Static', 'CFS', 'Static Noise', 'CFS Noise', 'location', 'northwest');
    legend boxoff
    title({subjects{i} '\pm1SEM (bootstrapped across runs)'});
    box off;
    pbaspect([1 1 1]);
    set(gcf, 'PaperPosition', [0 0 4 4])
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
    saveas(gcf, fullfile(outputDir, ['CFS_AssembleAMPlot_' subjects{i} '.pdf']), 'pdf');
    close(gcf);
end

%% Plot Group Average
averageAmplitudes = totalAmplitudes./length(subjects);
averageAmplitudeErrors = totalAmplitudeErrors./sqrt(length(subjects));
avgStaticBackground = averageAmplitudes(1);
avgCFSBackground = averageAmplitudes(2);
for j = 2:length(modulations)
    cfsData(j-1,:) = [averageAmplitudes((2*j-1):(2*j))]';
    cfsErrorData(j-1,:) = [averageAmplitudeErrors((2*j-1):(2*j))]';
end
graph = bar(gca,cfsData);
hold on;
set(gca,'XTick', 1:length(modulations)-1);
set(gca,'XTickLabel', modulations(2:4));
xlabel('Modulation');
set(gca, 'YLim', [0 6]);
    ylabel('Amplitude [\Delta%]');
title('Group Average');
plot(xlim, [avgStaticBackground avgStaticBackground], '--k');
plot(xlim, [avgCFSBackground avgCFSBackground], '--g');hold on;
set(graph,'BarWidth',1);    % The bars will now touch each other
numgroups = size(cfsData, 1);
numbars = size(cfsData, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));

for j = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, cfsData(:,j), cfsErrorData(:,j), 'k', 'linestyle', 'none');
end

legend('Static', 'CFS', 'Static Noise', 'CFS Noise', 'location', 'northwest');
legend boxoff

title({'Group average' '\pm1SEM (across subjects)'});

box off;
pbaspect([1 1 1]);
set(gcf, 'PaperPosition', [0 0 4 4])
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(outputDir, ['CFS_AssembleAMPlot_GrpAverage.pdf']), 'pdf');
close(gcf);
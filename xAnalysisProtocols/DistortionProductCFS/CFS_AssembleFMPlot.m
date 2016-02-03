function cfsDirectAttenuation
%% cfsDirectAttenuation

% Show amplitude response at both f and 2f
clear all; close all;
%% Set Subjects
subjects = {'s001' 's002' 's005'};
modulations = {'0.05 Hz' '0.5 Hz'};
outputDir = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions/AnalysisProtocols/DistortionProductCFS/Plots';
inputDir = '/Users/Shared/MATLAB/Experiments/OneLight/OLFlickerSensitivity/analysis/results/PupillometryCFSDirectedFM/';

hold all;
totalAmplitudes = zeros(length(modulations)*2+2,1);
totalAmplitudeErrors = zeros(length(modulations)*2+2,1); %used for computing group average
for j = 1:length(subjects)
    %% Open CSV and load subject data
    dataPath = strcat(inputDir, subjects(j),'-results.csv');
    data = csvread(char(dataPath),1,2);
    amplitudes = data(:,3)*100;
    amplitudeErrors = data(:,5)*100;
    totalAmplitudes = totalAmplitudes + amplitudes;
    totalAmplitudeErrors = totalAmplitudeErrors + amplitudeErrors;
    
    firstBackground = amplitudes(length(modulations)*2+1);
    secondBackground = amplitudes(length(modulations)*2+2);
    
    %% Plot Results
    for k = 1:length(modulations)
        cfsData(k,:) = [amplitudes(k) amplitudes(k + length(modulations))];
        cfsErrorData(k,:) = [amplitudeErrors(k) amplitudeErrors(k + length(modulations))];
    end
    graph = bar(cfsData, 0.2);
    hold on;
    set(gca,'XTick', 1:length(modulations));
    set(gca,'XTickLabel', modulations);
    set(gca, 'YLim', [0 6]);
    ylabel('Amplitude [\Delta%]');
    
    title(subjects(j));
    plot([0.6 1.4], [firstBackground firstBackground], '-k');
    plot([1.6 2.4], [secondBackground secondBackground], '-k');
    hold on;
    set(graph,'BarWidth',1);    % The bars will now touch each other
    numgroups = size(cfsData, 1);
    numbars = size(cfsData, 2);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for k = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, cfsData(:,k), cfsErrorData(:,k), 'k', 'linestyle', 'none');
    end
    
    title({subjects{j} '\pm1SEM (bootstrapped across runs)'});
    box off;
    pbaspect([1 1 1]);
    set(gcf, 'PaperPosition', [0 0 4 4])
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
    saveas(gcf, fullfile(outputDir, ['CFS_AssembleFMPlot_' subjects{j} '.pdf']), 'pdf');
    close(gcf);
    
end

legend('Static', 'CFS', 'location', 'southeast');

%% Plot Group Average
averageAmplitudes = totalAmplitudes./length(subjects);
averageAmplitudeErrors = totalAmplitudeErrors./sqrt(length(subjects));

avgFirstBackground = averageAmplitudes(length(modulations)*2+1);
avgSecondBackground = averageAmplitudes(length(modulations)*2+2);

for k = 1:length(modulations)
    cfsData(k,:) = [averageAmplitudes(k) averageAmplitudes(k + length(modulations))];
    cfsErrorData(k,:) = [averageAmplitudeErrors(k) ...
        averageAmplitudeErrors(k + length(modulations))];
end
graph = bar(cfsData, 0.2);
hold on;
set(gca,'XTick', 1:length(modulations));
set(gca,'XTickLabel', modulations);
set(gca, 'YLim', [0 6]);
ylabel('Increase in Amplitude (%)');
title('Group Average');
plot([0.6 1.4], [avgFirstBackground avgFirstBackground], '-k');
plot([1.6 2.4], [avgSecondBackground avgSecondBackground], '-k');
hold on;
set(graph,'BarWidth',1);    % The bars will now touch each other
numgroups = size(cfsData, 1);
numbars = size(cfsData, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));
for k = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, cfsData(:,k), cfsErrorData(:,k), 'k', 'linestyle', 'none');
end
title({'Group average' '\pm1SEM (across subjects)'});
box off;
pbaspect([1 1 1]);
set(gcf, 'PaperPosition', [0 0 4 4])
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(outputDir, ['CFS_AssembleFMPlot_GrpAverage.pdf']), 'pdf');
close(gcf);
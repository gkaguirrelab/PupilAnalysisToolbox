function PupilAnalysisToolbox_SummarizeDataQuality(Subjects, resultsPath)

averagesPath = '/Averages';

%% Create a director for this output, if it does not already exist
if ~exist(fullfile(resultsPath, averagesPath), 'dir');
    mkdir(fullfile(resultsPath, averagesPath));
end


%% Save the list of analyzed subjects
outFile = fullfile(resultsPath, averagesPath, 'subjectListforAverages.csv');
fid = fopen(outFile, 'w');
fprintf(fid, 'Subjects included in average plots:\n\n');
for s = 1:length(Subjects)
    fprintf(fid, [char(Subjects(s)),'\n']);
end
fclose(fid);

%% Load the vector of data quality from each subject

% suppress the warning string as tables are loaded
warningString='MATLAB:table:ModifiedVarnames';
warningS = warning('off',warningString);

for s = 1:length(Subjects)
    badTrialVectorOutFile = fullfile(resultsPath, [char(Subjects(s)) '/' char(Subjects(s)) '_PupilPulseData_BadTrialVector.csv']);
tmpTable=readtable(badTrialVectorOutFile);
trialProportionMissingData(s,:)=tmpTable.Prop_MissingData;
end

% restore the warning
warning(warningS.state,warningS.identifier);

%% Create the plot of across subject average data quality x trial
figure;
hold on;
shadedErrorBar(linspace(1,size(trialProportionMissingData,2),size(trialProportionMissingData,2)), ...
    mean(trialProportionMissingData, 1), ...
    std(trialProportionMissingData, 1)/sqrt(length(Subjects)));
pbaspect([1 1 1]);
xlabel('Trial number'); ylabel('Proportion missing data');
xlim([0 size(trialProportionMissingData,2)]); ylim([0 1]); 

title({['Mean\pm1SEM (n = ' num2str(length(Subjects)) ' subjects)']});

set(gca, 'TickDir', 'out');
set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 4 and height 4.
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 4 and height 4.
outFile1 = fullfile(resultsPath, averagesPath, ['DataQuality_x_Trial.pdf']);
saveas(gcf, outFile1, 'pdf');
outFile2 = fullfile(resultsPath, averagesPath, ['DataQuality_x_Trial.png']);
saveas(gcf, outFile2, 'png');
close(gcf);

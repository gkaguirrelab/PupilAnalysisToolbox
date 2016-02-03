basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/results/MelLightLevelDependence5sPulse_'


theDirs = {'MelanopsinDirectedPenumbralIgnoreNulledPositive' 'LMSDirectedNulledPositive'};
%% Load CSV files
theSubjects = {{'G092815Ax' 'J092915Rx' 'M092515Sx'} {'G100515Ax' 'J100615Rx' 'M100615Sx'} {'G092815Ax' 'J092915Rx' 'M092515Sx'} {'G100515Ax' 'J100615Rx' 'M100615Sx'} {'G092815Ax' 'J092915Rx' 'M092515Sx'}};

theNDvals = {'ND20' 'ND15' 'ND10' 'ND05' 'ND00'};
for n = 1:length(theNDvals)
    ndVal = theNDvals{n};
    for s = 1:length(theSubjects{n})
        tmpLMS{n, s} = csvread(fullfile([basePath ndVal], [char(theSubjects{n}{s}) ndVal '_PupilPulseData_LMS+_Mean.csv']));
        tmpMel{n, s} = csvread(fullfile([basePath ndVal], [char(theSubjects{n}{s}) ndVal '_PupilPulseData_Mel+_Mean.csv']));
        
        
    end
    
    set(gcf, 'PaperPosition', [0 0 8 4]); %Position plot at left hand corner with width 15 and height 6.
    set(gcf, 'PaperSize', [8 4]); %Set the paper to have width 15 and height 6.
    saveas(gcf, ['~/Desktop/MeanPupilSize_' ndVal '.pdf'], 'pdf');
    close(gcf);
end

%%
grandMeansLMS = cellfun(@mean,tmpLMS);
grandMeansMel = cellfun(@mean,tmpMel);
stdLMS = cellfun(@std,tmpLMS);
stdMel = cellfun(@std,tmpMel);
nLMS = cellfun(@length,tmpLMS);
nMel = cellfun(@length,tmpMel);
semLMS = stdLMS./sqrt(nLMS);
semMel = stdMel./sqrt(nMel);

offSet = 0.05;
theCols = [72 123 107 ; 238 159 96 ; 120 52 54 ]/255;
theLevels = 1:5;
for i = 1:3
    h = errorbar(theLevels-offSet, grandMeansLMS(:, i), semLMS(:, i), 'LineStyle', 'none', 'Color', 'k'); hold on;
    errorbar_tick(h, 0);
    plot(theLevels-offSet, grandMeansLMS(:, i), '-o', 'Color', theCols(i, :), 'MarkerFaceColor', theCols(i, :), 'MarkerEdgeColor', [0 0 0]);
    h = errorbar(theLevels+offSet, grandMeansMel(:, i), semMel(:, i), 'LineStyle', 'none', 'Color', 'k');
    errorbar_tick(h, 0)
    plot(theLevels+offSet, grandMeansMel(:, i), '-s', 'Color', theCols(i, :), 'MarkerFaceColor', theCols(i, :), 'MarkerEdgeColor', [0 0 0]);
    
end

pbaspect([1 1 1]);
ylim([2 7.5]); ylabel('Pupil diameter [mm]');

xlim([0 6]);
for i = 1:3
    h(i) = plot([1000 1001], [1000 1001], 's', 'Color', theCols(i, :), 'MarkerFaceColor', theCols(i, :), 'MarkerEdgeColor', [0 0 0])
end
h(4) = plot([1000 1001], [1000 1001], 'ok', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
h(5) = plot([1000 1001], [1000 1001], 'sk', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
legend(h, 'GA', 'JR', 'MS', 'L+M+S', 'Mel'); legend boxoff;

set(gca, 'XTick', 1:5, 'XTickLabel', {'ND2.0' 'ND1.5' 'ND1.0' 'ND0.5' 'ND0.0'}); xlabel('ND filter');
title({'Mean pupil size [5 sec before stimulus onset]', 'Error bars are \pm1SEM across trials'})


set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
saveas(gcf, '~/Desktop/MeanPupilSize.pdf', 'pdf');
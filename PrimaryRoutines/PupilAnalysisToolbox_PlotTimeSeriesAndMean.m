function PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath)

%%
ylimVal = 0.4;
for dd = 1:size(Data, 2)
    figure;
    hold on;
    plot([5 5], 100*[-ylimVal ylimVal], 'r');
    plot([min(Data(1, dd).t) max(Data(1, dd).t)], [0 0], '-k');
    shadedErrorBar(Data(1, dd).t, 100*nanmean([Data(:, dd).AvgTimeSeries], 2), 100*nanstd([Data(:, dd).AvgTimeSeries], [], 2)/sqrt(size([Data(:, dd).AvgTimeSeries], 2)));
    xlabel('Time [s]');
    ylabel('\Delta%');
    pbaspect([1 1 1]);
    title({strrep(Data(1, dd).label, '_', ' ')});
    xlim([min(Data(1, dd).t) max(Data(1, dd).t)]);
    ylim([-ylimVal ylimVal]);
    set(gca, 'TickDir', 'out');
    set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
    saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(1, dd).label, '_', ' ') '.pdf']), 'pdf');
    close(gcf);
end
%%
figure;
subjectCols = [224 236 244 ; 158 188 218 ; 136 86 167]/255;
nSubjects = size(Data, 1);
nConds = size(Data, 2);
shiftVals = [-0.2 0 0.2];
for dd = 1:size(Data, 2)
    hold on;
    for ii = 1:nSubjects
       tmp = plot(Data(ii, dd).Mean, dd+shiftVals(ii)+rand(size(Data(ii, dd).Mean))/50, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', subjectCols(ii, :), 'MarkerSize', 8);
       h(ii) = tmp(1);
       plot([mean(Data(ii, dd).Mean) mean(Data(ii, dd).Mean)], [dd+shiftVals(ii)+0.5*shiftVals(1) dd+shiftVals(ii)+0.5*shiftVals(end)], '-', 'Color', 'r', 'LineWidth', 3);
    end
    xlabel('Diameter [mm]');
    ylabel('Condition');
    pbaspect([1 1 1]);
    xlim([1 9]); ylim([0.5 nSubjects+2.5]);
    title('Pupil diameter mean');
end
set(gca, 'YTick', 1:nConds, 'YTickLabel', {Data(1, :).label});
legend(h, Subjects, 'Location', 'NorthEast'); legend boxoff;

set(gca, 'TickDir', 'out');
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 15 and height 6.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 15 and height 6.
saveas(gcf, fullfile(resultsPath, 'MeanPupilDiameter.pdf'), 'pdf');
close(gcf);

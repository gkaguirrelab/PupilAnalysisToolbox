function PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params)
%  PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params)

for dd = 1:size(Data, 2)
    figure;
    hold on;
    plot([0 params.xLim], [0 0], '-', 'Color', [0.3 0.3 0.3]); hold on;
    shadedErrorBar(Data(1, dd).timeSecs, nanmean([Data(:, dd).AvgTimeSeries], 2), ...
        nanstd([Data(:, dd).AvgTimeSeries], [], 2)/sqrt(size([Data(:, dd).AvgTimeSeries], 2)));
    plot([0 params.PulseDurationSecs], [0.1 0.1], '-r', 'LineWidth', 2);
    pbaspect([1 1 1]);
    xlim([-params.PulseOnsetSecs params.xLim]); ylim([-params.yLim params.yLim]);
    xlabel('Time [Secs]'); ylabel('Pupil diameter [change]');
    title({strrep(Data(1, dd).label, '_', ' ') ['Mean\pm1SEM (n = ' num2str(size(Data, 1)) ' subjects)']});
    
    set(gca, 'TickDir', 'out');
    set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 4 and height 4.
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 4 and height 4.
    outFile1 = fullfile(resultsPath, ['TimeSeries_' strrep(Data(1, dd).label, '_', ' ') '.pdf']);
    saveas(gcf, outFile1, 'pdf');
    outFile2 = fullfile(resultsPath, ['TimeSeries_' strrep(Data(1, dd).label, '_', ' ') '.png']);
    saveas(gcf, outFile2, 'png');
    close(gcf);
    
end

%%
% for jj = 1:size(Data, 1)
%     diff(:, jj) = Data(jj, 2).AvgTimeSeries-Data(jj, 3).AvgTimeSeries;
% end
% hold on
%     plot([5 5], 100*[-params.yLim params.yLim], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 2).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 2).label, '_', ' ') ' - ' strrep(Data(1, 3).label, '_', ' ')]);
%     xlim([0 params.xLim]);
%     ylim(100*[-params.yLim params.yLim]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, 2).label, '_', ' ') ' - ' strrep(Data(3, 3).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);

% for jj = 1:size(Data, 1)
%     diff(:, jj) = Data(jj, 2).AvgTimeSeries-Data(jj, 1).AvgTimeSeries;
% end
% hold on
%     plot([5 5], 100*[-params.yLim params.yLim], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 2).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 2).label, '_', ' ') ' - ' strrep(Data(1, 1).label, '_', ' ')]);
%     xlim([0 params.xLim]);
%     ylim(100*[-params.yLim params.yLim]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, 2).label, '_', ' ') ' - ' strrep(Data(3, 1).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);
%
%     for jj = 1:size(Data, 1)
%     diff(:, jj) = Data(jj, 3).AvgTimeSeries-Data(jj, 1).AvgTimeSeries;
% end
% hold on
%     plot([5 5], 100*[-params.yLim params.yLim], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 3).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 3).label, '_', ' ') ' - ' strrep(Data(1, 1).label, '_', ' ')]);
%     xlim([0 params.xLim]);
%     ylim(100*[-params.yLim params.yLim]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, 3).label, '_', ' ') ' - ' strrep(Data(3, 1).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);
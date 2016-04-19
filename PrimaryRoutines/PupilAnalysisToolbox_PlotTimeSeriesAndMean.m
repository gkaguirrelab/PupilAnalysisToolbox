function PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, xLimVal, yLimVal)

%%
% for dd = 1:size(Data, 2)
%     figure;
%     hold on;
%     plot([5 5], 100*[-yLimVal yLimVal], 'r');
%     plot([min(Data(1, dd).t) max(Data(1, dd).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, dd).t, 100*nanmean([Data(:, dd).AvgTimeSeries], 2), 100*nanstd([Data(:, dd).AvgTimeSeries], [], 2)/sqrt(size([Data(:, dd).AvgTimeSeries], 2)));
%     xlabel('Time [s]');
%     ylabel('\Delta%');
%     pbaspect([1 1 1]);
%     title({strrep(Data(1, dd).label, '_', ' ')});
%     xlim([0 xLimVal]);
%     ylim(100*[-yLimVal yLimVal]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(1, dd).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);
% end

%%
% for jj = 1:size(Data, 1)
%     diff(:, jj) = Data(jj, 2).AvgTimeSeries-Data(jj, 3).AvgTimeSeries;
% end
% hold on
%     plot([5 5], 100*[-yLimVal yLimVal], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 2).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 2).label, '_', ' ') ' - ' strrep(Data(1, 3).label, '_', ' ')]);
%     xlim([0 xLimVal]);
%     ylim(100*[-yLimVal yLimVal]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, 2).label, '_', ' ') ' - ' strrep(Data(3, 3).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);

% for jj = 1:size(Data, 1)
%     diff(:, jj) = Data(jj, 2).AvgTimeSeries-Data(jj, 1).AvgTimeSeries;
% end
% hold on
%     plot([5 5], 100*[-yLimVal yLimVal], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 2).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 2).label, '_', ' ') ' - ' strrep(Data(1, 1).label, '_', ' ')]);
%     xlim([0 xLimVal]);
%     ylim(100*[-yLimVal yLimVal]);
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
%     plot([5 5], 100*[-yLimVal yLimVal], 'r');
%     plot([min(Data(1, 1).t) max(Data(1, 1).t)], [0 0], '-k');
%     shadedErrorBar(Data(1, 3).t, 100*nanmean(diff, 2), 100*nanstd(diff, [], 2)/sqrt(size(diff, 2)));
%     xlabel('Time [s]');
%     ylabel('Difference in \Delta%');
%     pbaspect([1 1 1]);
%     title([strrep(Data(1, 3).label, '_', ' ') ' - ' strrep(Data(1, 1).label, '_', ' ')]);
%     xlim([0 xLimVal]);
%     ylim(100*[-yLimVal yLimVal]);
%     set(gca, 'TickDir', 'out');
%     set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
%     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
%     saveas(gcf, fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, 3).label, '_', ' ') ' - ' strrep(Data(3, 1).label, '_', ' ') '.pdf']), 'pdf');
%     close(gcf);
function PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, comparisonsIdx, comparisonsPolarity)
% PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params)

%% Plot the averages
for dd = 1:size(Data, 2)
    figure;
    hold on;
    plot([-params.PulseOnsetSecs params.xLim], [0 0], '-', 'Color', [0.3 0.3 0.3]); hold on;
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

%% Plot comparisons
for ii = 1:length(comparisonsIdx)
    for jj = 1:size(Data, 1)
        diff(:, jj) = Data(jj, comparisonsPolarity{ii}(1)*comparisonsIdx{ii}(1)).AvgTimeSeries + comparisonsPolarity{ii}(2)*Data(jj, comparisonsIdx{ii}(2)).AvgTimeSeries;
    end
    hold on
    plot([-params.PulseOnsetSecs params.xLim], [0 0], '-', 'Color', [0.3 0.3 0.3]); hold on;
    shadedErrorBar(Data(1, 2).timeSecs, nanmean(diff, 2), nanstd(diff, [], 2)/sqrt(size(diff, 2)));
    plot([0 params.PulseDurationSecs], [0.1 0.1], '-r', 'LineWidth', 2);
    pbaspect([1 1 1]);
    xlim([-params.PulseOnsetSecs params.xLim]); ylim([-params.yLim params.yLim]);
    xlabel('Time [Secs]'); ylabel('Pupil diameter [change]');
    title({[strrep(Data(2, comparisonsIdx{ii}(1)).label, '_', ' ') '-' strrep(Data(3, comparisonsIdx{ii}(2)).label, '_', ' ')] ['Mean\pm1SEM (n = ' num2str(size(Data, 1)) ' subjects)']});
    set(gca, 'TickDir', 'out');
    set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
    set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
    outFile1 = fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, comparisonsIdx{ii}(1)).label, '_', ' ') '-' strrep(Data(3, comparisonsIdx{ii}(2)).label, '_', ' ') '.pdf']);
    outFile2 = fullfile(resultsPath, ['TimeSeries_' strrep(Data(2, comparisonsIdx{ii}(1)).label, '_', ' ') '-' strrep(Data(3, comparisonsIdx{ii}(2)).label, '_', ' ') '.png']);
    saveas(gcf, outFile1, 'pdf');
    saveas(gcf, outFile2, 'png');
    close(gcf);
end
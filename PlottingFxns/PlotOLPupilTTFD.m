function PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,Data,ErrorVals,PlotTitle,yLabelText,yRange,figureID,subplotIndex,maxSubplots)
% Plot the TTFD for all directions

% If the JITTER flag is set the 1, points on x axis are jittered to make
% the plot less cluttered
JITTER = 0;
if JITTER
    jitterFactor = linspace(-0.05, 0.05, length(UniqueDirectionLabels));
else
    jitterFactor = zeros(1, length(UniqueDirectionLabels));
end

[ Colors, MarkerStyles, LineStyles ] = PlotDirectionColors( UniqueDirectionLabels );

figure(figureID);


if (maxSubplots == 1)
    subplot(1, 1, subplotIndex);
end

if (maxSubplots == 2)
    subplot(2, 1, subplotIndex);
end

if (maxSubplots == 3)
    subplot(3, 1, subplotIndex);
end

if (maxSubplots == 4)
    subplot(2, 2, subplotIndex);
end

if (maxSubplots == 5)
    subplot(2, 3, subplotIndex);
end

if (maxSubplots == 6)
    subplot(2, 3, subplotIndex);
end


hold on;
if (not(isempty(PlotTitle)))
    title(char(PlotTitle));
end
if (not(isempty(yRange)))
    ylim(yRange);
end


xlim([-2.5 .5]);
xlabel(['log_{10} frequency [Hz]']);
ylabel(yLabelText);
if (isempty(ErrorVals))
    for d = 1:length(UniqueDirectionLabels)
        plot([log10(UniqueFreqs)]+jitterFactor(d),Data(:,d),'Color', Colors(d, :), 'Marker', MarkerStyles{d}, 'LineStyle', LineStyles{d},'MarkerSize',5, 'MarkerFaceColor', Colors(d, :), 'MarkerEdgeColor', Colors(d, :));
    end
else
    for d = 1:length(UniqueDirectionLabels)
        ph = errorbar([log10(UniqueFreqs)]+jitterFactor(d),Data(:,d),ErrorVals(:,d),'Color', Colors(d, :), 'Marker', MarkerStyles{d}, 'LineStyle', LineStyles{d},'MarkerSize',5, 'MarkerFaceColor', Colors(d, :), 'MarkerEdgeColor', Colors(d, :));
    end
end



% Plot 0 line
plot([-2.5 .5], [0 0], '--k')

if maxSubplots == subplotIndex
    legend(UniqueDirectionLabels);
end
legend boxoff;

pbaspect([1 1 1]);

end
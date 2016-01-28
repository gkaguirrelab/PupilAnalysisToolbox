function figureID = PlotOLPupilSparklines( DataMatrixIn, FitMatrixIn, ErrorMatrixIn, ErrorMultiplier, yRange, figureID_in, RowLabels, ColumnPlotTitle, ColumnIndex, NumColumns )


ErrorColor=[0.6    0.6    0.6]; % light gray
FitColor=[1    0.0000    0.000];  % red
DataColor=[.3    .3    .3]; % dark gray

DataMatrix=DataMatrixIn;
FitMatrix=FitMatrixIn;
ErrorMatrix=ErrorMatrixIn;

tmpSize=size(DataMatrix);
NumSeries=tmpSize(1);
DataLength=tmpSize(2);

% A kludge to detect when the rows and columns of the passed data have been
% swaped, and to correct this.

if (NumSeries > DataLength)
    DataMatrix=DataMatrix';
    FitMatrix=FitMatrix';
    ErrorMatrix=ErrorMatrix';
    tmp=DataLength;
    DataLength=NumSeries;
    NumSeries=tmp;
end

if ~isempty(ErrorMatrixIn)
    PlusErrorMatrix=DataMatrix + (ErrorMultiplier*ErrorMatrix);
    NegErrorMatrix=DataMatrix - (ErrorMultiplier*ErrorMatrix);
else
    PlusErrorMatrix=DataMatrix;
    NegErrorMatrix=DataMatrix ;
end

if (isempty(figureID_in))
    figureID=figure('Color',[1 1 1]);
else
    figureID=figureID_in;
end

if not(isempty(ColumnPlotTitle))
    ColumnWidth=1/NumColumns;
    Position=[ColumnWidth*(ColumnIndex-1) 0.8 ColumnWidth 0.1];
    annotation('textbox', Position, ...
        'String', allwords(ColumnPlotTitle, '_'), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',10)
end % if placing title

for ni=1:NumSeries
    subplot_tight(NumSeries+1,NumColumns, (ni)*NumColumns+ColumnIndex);
    hold on;
    ylim(yRange);
    plot([1 DataLength], [0 0], '--b')              % plot the zero line
    if and((ni==1),(ColumnIndex==1))
        plot([0 0],yRange,'-b') % add a scale bar to the first column
    end
    
    % Only plot the vector if it is not entirely composed of zeros
    
    if ~(and(max(PlusErrorMatrix(ni,:))==0,min(PlusErrorMatrix(ni,:))==0))
        plot(PlusErrorMatrix(ni,:),'Color', ErrorColor, 'LineStyle', '-');
        plot(NegErrorMatrix(ni,:),'Color', ErrorColor, 'LineStyle', '-');
    end
    
    if ~(and(max(FitMatrix(ni,:))==0,min(FitMatrix(ni,:))==0))
        plot(FitMatrix(ni,:),'Color', FitColor, 'LineStyle', '-');
    end
    
    if ~(and(max(DataMatrix(ni,:))==0,min(DataMatrix(ni,:))==0))
        plot(DataMatrix(ni,:),'Color', DataColor, 'LineStyle', '-');
    end
    
    set(gca, 'Visible', 'off');                     % Turn off the axes
    set(gca, 'Color', [0.1 0.1 0.1]);                     % Turn off the axes
end



hold off

end


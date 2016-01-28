function [] = PlotOLPupilPolar(Amplitudes,Phases,AmpError,PhaseError,MaxAmp,PlotTitle,figureID,subplotIndex,maxSubplots,UniqueDirectionLabels,varargin)

if ~isempty(varargin)
   hugeMarker = varargin{1}; 
end

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

n_data = length(Amplitudes);

fake = polar(Phases,MaxAmp*ones(1,n_data)); set(fake,'Visible','off'); hold on;
if (not(isempty(PlotTitle)))
    title(char(PlotTitle));
end

th = findall(gcf,'Type','text');

for i = 1:length(th),
    set(th(i),'FontSize',12)
end

[ Colors, MarkerStyles, LineStyles ] = PlotDirectionColors( UniqueDirectionLabels );

for ni = n_data:-1:1
    h=polar(Phases(ni),Amplitudes(ni),char(MarkerStyles(ni))); hold on;
    if ~isempty(varargin)
    set(h,'markersize',12);
    else
        set(h,'markersize',5);
    end
       
    set(h,'markerfacecolor',Colors(ni,:));
    set(h,'LineWidth',1);
    set(h,'markeredgecolor',Colors(ni,:));
    
    % Add error bars in the amplitude and then phase directions
    
    if (Amplitudes(ni)-AmpError(ni)) < 0
        h=polar(Phases(ni)*ones(1,3),[0, Amplitudes(ni), Amplitudes(ni)+AmpError(ni)],'-r'); hold on;
        set(h,'LineWidth',.75);
    else
        h=polar(Phases(ni)*ones(1,3),[Amplitudes(ni)-AmpError(ni), Amplitudes(ni), Amplitudes(ni)+AmpError(ni)],'-r');  hold on;
        set(h,'LineWidth',.75);
    end
    
    ErrorBarPhases=linspace((Phases(ni)-PhaseError(ni))*10,(Phases(ni)+PhaseError(ni))*10,round(PhaseError(ni)*20))/10;
    ErrorBarAmps=Amplitudes(ni)*ones(1,length(ErrorBarPhases));
    
    h=polar( ErrorBarPhases ,ErrorBarAmps,'-r');  hold on;
    set(h,'LineWidth',.75);
    
end

hold off

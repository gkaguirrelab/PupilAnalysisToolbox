function [ Colors, MarkerStyles, LineStyles ] = PlotDirectionColors( UniqueDirectionLabels )

%PLOTDIRECTIONCOLORS Standardizes plot colors for different modulation
%  directions

NumDirections=length(UniqueDirectionLabels);

Colors=ones(NumDirections,3)*(-1);


% Implement a look-up table for direction labels using an embarassing set
% of if statements

for i=1:NumDirections
    if strcmp(UniqueDirectionLabels(i),'Iso')
        Colors(i,:)=[0.4784    0.1843    0.3020]; % brown
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    if strcmp(UniqueDirectionLabels(i),'LM')
        Colors(i,:)=[0.9490    0.7216    0.0275]; % yellow
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LMDirected02')
        Colors(i,:)=[0.3 0 0]; % red
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LMDirected04')
        Colors(i,:)=[0.5 0 0]; % red
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LMDirected08')
        Colors(i,:)=[0.7 0 0]; % red
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LMDirected32')
        Colors(i,:)=[1 0 0]; % red
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
                    
    if strcmp(UniqueDirectionLabels(i),'Mel')
        Colors(i,:)=[0.0275    0.4902    0.9490]; % cyan
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    
    if strcmp(UniqueDirectionLabels(i),'MelanopsinDirected')
        Colors(i,:)=[0.0275    0.4902    0.9490]; % cyan
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'MelScaled')
        Colors(i,:)=[0.1275    0.5902    0.8390]; % cyan
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'S')
        Colors(i,:)=[0    0.0588    0.5098];  % blue
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'SDirected')
        Colors(i,:)=[0    0.0588    0.5098];  % blue
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    
    if strcmp(UniqueDirectionLabels(i),'S_r')
        Colors(i,:)=[0    0.0588    0.5098];  % blue
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    if strcmp(UniqueDirectionLabels(i),'L/M + Mel')
        Colors(i,:)=[0.8667    1    0.0667];  % Neon
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'L/M + Mel + S')
        Colors(i,:)=[0    0    0];  % black
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    if strcmp(UniqueDirectionLabels(i),'Mel + S')
        Colors(i,:)=[0    0    0];  % black
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    if strcmp(UniqueDirectionLabels(i),'L/M + S')
        Colors(i,:)=[0    0    0];  % black
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    
    
    if strcmp(UniqueDirectionLabels(i),'S/Mel-IN')
        Colors(i,:)=[0.4784    0.0627    0.8941];  % purple
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'S/Mel-OUT')
        Colors(i,:)=[0.4784    0.0627    0.8941];  % purple
        MarkerStyles(i)={'p'}; % pentagram
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'Noise')
        Colors(i,:)=[0.5    0.5    0.5];  % gray
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LM-Decrease')
        Colors(i,:)=[0.9490    0.7216    0.0275]; % yellow
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'LM-Increase')
        Colors(i,:)=[0.9490    0.7216    0.0275]; % yellow
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'S-Decrease')
        Colors(i,:)=[0    0.0588    0.5098];  % blue
        MarkerStyles(i)={'s'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'S-Increase')
        Colors(i,:)=[0    0.0588    0.5098];  % blue
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    
    if strcmp(UniqueDirectionLabels(i),'RodDirected')
        Colors(i,:)=[0.5451    0.4667    0.3961];  % blue
        MarkerStyles(i)={'o'};
        LineStyles(i)={'-'};
    end
    
    if strcmp(UniqueDirectionLabels(i),'OmniSilent')
        Colors(i,:)=[0.5 1 0];  % green
        MarkerStyles(i)={'^'};
        LineStyles(i)={'-'};
    end
    
    if (Colors(i,1)==-1)
        Colors(i,:)=[0    0    0];  % black
        MarkerStyles(i)={'x'};
        LineStyles(i)={'-'};
    end % catch the condition of no matches to labels
    
end % embarrasing loop of case construction


end


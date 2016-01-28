function [AvgPS]=MakeAvgPS(TimeSeries,TrialFrequencies,TrialDirections,UniqueDirectionLabels,sampling_frequency,SubjectLabel)

% This routine implements a control analysis for the TTF4D data for the
% PNAS paper. The goal is to measure the average power spectrum of the
% residual pupil data across trials for each subject. We cannot use
% tradditional FFT methods to calculate the power spectrum, as we have
% missing data (NaNs) in our time series. Therefore, we adopt the
% approach of using the same routine we use to fit our data and measuring
% the amplitude of variation in pupil size at the set of frequencies
% between 0.01 Hz to HighestFreq (set to be 2 Hz). This is done for
% all frequencies in a trial except those frequencies that are an odd
% or even harmonic of the stimulation frequency, thus providing the
% power spectrum of the residual.
%
% A wrinkle is that in the TTF4D data, the stimulus was not a perfect
% sinusoid.  The imperfection of the stimulus sinusoids create side-lobes
% around the targeted frequency. To handle this, we exclude from
% the residual power spectrum not only the harmonics, but the frequencies
% adjacent to the harmonics.
%
% - GKA, 9/8/2014



fprintf('  - Calculating average power spectrum.\n'); % Notify user

HighestFreq=2;  % Highest frequency in the plot, expressed in Hz.

sizer=size(TimeSeries);
PSMatrix=zeros(sizer(1),200);

for trial=1:sizer(1)
    x=squeeze(TimeSeries(trial,:));
    for freq=1:200
        IsStimFreq=0;
        
        f=freq/(2000/sampling_frequency);
        if (mod(f,TrialFrequencies(trial)))==0
            IsStimFreq=1;
        end
        
        fplus=(freq+1)/(2000/sampling_frequency);
        if (mod(fplus,TrialFrequencies(trial)))==0
            IsStimFreq=1;
        end
        
        if freq ~= 1
            fminus=(freq-1)/(2000/sampling_frequency);
            if (mod(fminus,TrialFrequencies(trial)))==0
                IsStimFreq=1;
            end
        end
        
        if IsStimFreq==1
            Amp=NaN;
        else
            [Amp,~,~] =...
                LeastSquaresSpectralFit(x,f,sampling_frequency);
        end
        PSMatrix(trial,freq)=Amp;
    end
end

for freq=1:200
    f(freq)=freq/(2000/sampling_frequency);
end

AvgPS=nanmean(PSMatrix,1);
NANSTD_NORMFLAG=1;
StdvPS=nanstd(PSMatrix,NANSTD_NORMFLAG,1);

figAvgPS = figure;
plot(log10(f),100*(AvgPS+StdvPS),'-r');
hold on
plot(log10(f),100*(AvgPS),'-k');
plot(log10(f),100*(AvgPS-StdvPS),'-r');

plotcolors=['-m','-y','-b','-c'];

directions=unique(TrialDirections);

for d=1:length(directions)
    Indices = find(strcmp(TrialDirections,directions(d)));
    
    % Create the average time series. This is done by first
    %  assembling a matrix across time-series, and then using
    %  the nanmean to handle NaN points
    
    for i=1:length(Indices)
        if (i==1)
            SubMatrix=PSMatrix(Indices(i),:);
        else
            SubMatrix(i,:)=PSMatrix(Indices(i),:);
        end % if the first timeseries
    end % for number of indicies
    
    subAvgPS=nanmean(SubMatrix,1);
    
    plot(log10(f),100*subAvgPS,plotcolors(d));
end

xlabel(['log_{10} frequency [Hz]']);
ylabel(['mean aross trial absolute pupil amplitude [% change]']);
PlotTitle=['Subject ' char(SubjectLabel) ' - Avg PS of residuals [±1SD] '];
title(char(PlotTitle));
saveas(figAvgPS, [char(SubjectLabel) '-AvgPS'], 'pdf');

hold off

fprintf('  - Done with average power spectrum.\n'); % Notify user


end
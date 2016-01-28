function [NoiseMeanAmplitude,NoiseMeanPhase,NoiseSEMAmplitude,NoiseSEMPhase]=NoiseBootStrap(TimeSeries,IndicesToSample,TestFrequency,sampling_frequency)

% TimeSeries -- The array of cleaned time-series data from each trial
% IndicesToSample -- An array of array indices in TimeSeries that
%   correspond to the particular temporal frequency and modulation direction
%   which is being assessed
% TestFrequency -- The temporal frequeuncy of the stimulus modulation
% sampling_frequency -- the sampling frequency of the TimeSeries data.

% This routine calculates the standard deviation of bootstraps where each
% is the average of the full n trials of data, sampled with a random sample
% with replacement. This approximates the standard error of the mean.


% Test for Matlab version. The routine requires "randperm" which is not
% available in early versions. If Matlab is too old, just return zeros.

if (strcmp(version,'7.11.0.584 (R2010b)'))
    NoiseMeanAmplitude=0;
    NoiseMeanPhase=0;
    NoiseSEMAmplitude=0;
    NoiseSEMPhase=0;
    fprintf('Matlab version too old to implement bootstrap error calculation');
else  % proceed with bootstrap
    
    
    % VARIABLE DECLARATION
    
    nBootstraps = 1000;
    
    lengthIndicesToSample = length(IndicesToSample);
    
    BootstrapAmplitudes = zeros(1,nBootstraps);
    BootstrapPhases = zeros(1,nBootstraps);
    
    for x = 1:nBootstraps
        NoiseTrialAverage=zeros(1,2000);
        
        % This will sample with replacement from the indices of the trials,
        % matching the total number of available trials.
        
        SampleWithReplace= ...
            datasample(linspace(1,lengthIndicesToSample, ...
            lengthIndicesToSample),lengthIndicesToSample,'Replace',true);
        
        RandomSampleIndices = IndicesToSample(SampleWithReplace);
        
%       Create the average time series. This is done by first
%       assembling a matrix across time-series, and then using
%       the nanmean to handle NaN points
                
        for i=1:lengthIndicesToSample
            if (i==1)
                TimeSeriesMatrix=[TimeSeries(RandomSampleIndices(i),:)'];
            else
                TimeSeriesMatrix=[TimeSeriesMatrix TimeSeries(RandomSampleIndices(i),:)'];
            end % if the first timeseries
        end % for number of samples
        NoiseTrialAverage=nanmean(TimeSeriesMatrix,2)';
        
        [BootstrapAmplitudes(x),BootstrapPhases(x),BestFit_Discard] = LeastSquaresSpectralFit(NoiseTrialAverage,TestFrequency,sampling_frequency);
    end % across Bootstraps
    
    % Note that when we take the standard deviation across the bootstraps
    % we are obtaining an estimate of the standard error of the mean in the
    % full length data. 
    
    NoiseMeanAmplitude=mean(BootstrapAmplitudes);
    NoiseSEMAmplitude=std(BootstrapAmplitudes);
    [NoiseMeanPhase,~,~]=circ_mean(BootstrapPhases,[],[]);
    [~,NoiseSEMPhase]=circ_std(BootstrapPhases,[],[],[]);
    
end % if statement checking for Matlab version

end

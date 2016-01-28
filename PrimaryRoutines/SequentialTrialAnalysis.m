function [TimeSeries]=SequentialTrialAnalysis(...
    sampling_frequency,...
    full_trial_length,...
    final_trial_length,...
    adapt_length,...
    minimum_length_trial,...
    StutterErrorFlag,...
    spike_remover_params,...
    sgolay_span,...
    sgolay_polynomial,...
    BadPercentChangeThreshold,...
    BadNanThreshold,...
    Subjects,...
    Protocols,...
    SyntheticIsoFlag,...
    SaveDataFlag,...
    SavePlotsFlag,...
    RelabelDirections,...
    newLabels,...
    oldLabels,...
    ResultsDirName,...
    configFileNames,...
    whichSttingIndexToValidate,...
    TrialInspectorFlag,...
    HarmonicModelFlag,...
    HarmonicTestFlag,...
    StimOnsetDelay,...
    RelabelSubjectsFlag,...
    StackPlotFlag,...
    OPN4ScaleFlag,...
    OLFlickerSensitivityFlag,...
    analyzeGazeData)


% Set a default value for max. contrast. We only use this for CRF
% measurements.
maxContrast = 0.5;


if ~exist('analyzeGazeData', 'var')
    analyzeGazeData = false;
end


if exist('OLFlickerSensitivity', 'var');
    OLFlickerSensitivity = false;
end

% Add the OLPupilDiameter folder to the path. Different conventions.
if isdir('/Users/Shared/Matlab/Experiments/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OLPupilDiameter';
elseif isdir('/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter';
end

% Add the OLPupilDiameter folder to the path. Different conventions.
if isdir('/Users/Shared/Matlab/Experiments/OLFlickerSensitivity/');
    basePath = '/Users/Shared/Matlab/Experiments/OLFlickerSensitivity';
elseif isdir('/Users/Shared/Matlab/Experiments/OneLight/OLFlickerSensitivity/');
    basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLFlickerSensitivity';
end


ResultsFullPath = fullfile(basePath, 'analysis', 'results');

% Calculate the  length of our data vectors
datalen = full_trial_length*sampling_frequency;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RETRIEVE CONTRASTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to obtain the contrasts that each of the photoreceptors gets.
% This is done with a call to OLPDCalculateContrast(...).
configPath = fullfile(basePath, 'code', 'config');

% Load the config files. We hard code 0 in there. This is the index for
% which we want to calculate contrasts. That is 90 deg phase angle when we
% have a stimulus titration with 200 steps between 0 and 2*pi.
%

if ~(isempty(configFileNames))
    for i = 1:length(configFileNames)
        [contrastTarget(i, :) contrastPredicted(i, :) receptorNames] = OLPDCalculateContrasts(configFileNames{i}, false, whichSttingIndexToValidate, false);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ROUTINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now loop across subjects and perform the analysis

for SubjectID=1:length(Subjects)
    
    % Close all existing subject-specific plots
    close all;
    
    % These figures will contain results for each subject
    figMeanPupilByDirection = figure;
    figNoiseByDirection = figure;
    figPhase = figure;
    figAmplitude = figure;
    figBackgroundAdapt = figure;
    
    if (SyntheticIsoFlag==1)
        figComboScaledPolarPlots = figure;
    end
    
    if (TrialInspectorFlag==1)
        figTrialInspector=figure;
        pause on;
    end
    
    figPolarPlots = figure;
    figSparkLinePlots = figure;
    figCycleAveragePlots = figure;
    if (SyntheticIsoFlag==1)
        figScaledPolarPlots = figure;
    end
    
    % Clear out variables from the prior subject
    clear Data;
    clear TimeSeries;
    clear TimeSeriesMatrix;
    clear TrialDirections;
    clear TrialFrequencies;
    clear TrialPhases;
    clear AmplitudesByTrial;
    clear AvgModelFits;
    clear AvgModelFits;
    clear BackgroundAdaptData;
    clear TrialDateTime;
    
    % Clear the figure windows used for each subject
    
    if (TrialInspectorFlag==1)
        figure(figTrialInspector);
        clf;
    end
    
    figure(figPolarPlots);
    clf;
    figure(figSparkLinePlots);
    clf;
    figure(figCycleAveragePlots);
    clf;
    if (SyntheticIsoFlag==1)
        figure(figScaledPolarPlots);
        clf;
    end
    
    if (RelabelSubjectsFlag==1)
        SubjectLabel = ['S' num2str(SubjectID)];
    else
        SubjectLabel=Subjects(SubjectID);
    end
    
    fprintf(['\n> Loading and concatenating sessions for subject '...
        char(Subjects(SubjectID))]); % Notify user
    
    % Loop over sessions and directions for a subject and concatenate the data
    
    FirstGoodSessionFlag=1;
    FirstStageBadTrials=0;
    
    for sess=1:10
        for p=1:length(Protocols)
            
            % Load data set
            
            if OLFlickerSensitivityFlag
                [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempAlphaSpacing,TempDateTime]=...
                    DataLoaderFlickerSensitivity(char(Protocols(p)),char(Subjects(SubjectID)),sess);
            else
                [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempAlphaSpacing,TempDateTime]=...
                    DataLoader(char(Protocols(p)),char(Subjects(SubjectID)),sess);
            end
            
            % Reshape the data if we have to. If data collection is aborted
            % and saved out manually, the data struct will likely be wrong.
            % We do a simple check to see if the dimensions are right (it
            % should be (Nx1), and transpose the matrix if not.
            if size(TempData, 1) < size(TempData, 2)
                TempData = TempData';
            end
            
            % Iterate over TempData to fill in empty diameter/time fields
            % with NaN
            for a = 1:length(TempData)
               if isempty(TempData(a).diameter)
                   TempData(a).diameter = NaN;
                   TempData(a).time = NaN;
               end
            end
           
            
            TempTrialFirstSample = cell2mat(arrayfun(@(s) s.time(1), TempData,'uni',0))';
            %if isempty(TempData.time_inter)
            %    TempTrialFirstSampleInter = 0;
            %else
            TempTrialFirstSampleInter = cell2mat(arrayfun(@(s) s.time_inter(1), TempData,'uni',0))';
            
            % Sometimes all of our data were good. Then, the time stamp of
            % the first interruption remained at 0. That's <960-980, so we
            % set it to Inf.
            TempTrialFirstSampleInter(TempTrialFirstSampleInter == 0) = Inf;
            TempTrialFirstSample = min([TempTrialFirstSample ; TempTrialFirstSampleInter]);
            TempDateTime = repmat({TempDateTime}, 1, length(TempData));
            
            if (not(isempty(TempData)))
                
                fprintf('.'); % Notify user
                
                % Keep AlphaSpacing
                
                AlphaSpacing=TempAlphaSpacing;
                
                % The first two trials are discarded from the primary analysis for each session, as they
                % correspond to the adaptation period and the initial wrap-around trial
                % The first trial is kept, however, and averaged for the
                % subject across all directions to provide the response for
                % adaptation to the background.
                
                if (FirstGoodSessionFlag==1)
                    BackgroundAdaptData=TempData(1);
                else
                    BackgroundAdaptData=[BackgroundAdaptData TempData(1)];
                end
                
                TempData(2)=[]; TempTrialFrequencies(2)=[]; TempTrialDirections(2)=[]; TempTrialPhases(2)=[]; TempTrialFirstSample(2) = []; TempTrialFirstSampleInter(2) = []; TempDateTime{2} = []; TempDateTime(cellfun(@(s) isempty(s),TempDateTime))=[];
                TempData(1)=[]; TempTrialFrequencies(1)=[]; TempTrialDirections(1)=[]; TempTrialPhases(1)=[]; TempTrialFirstSample(1) = []; TempTrialFirstSampleInter(1) = []; TempDateTime{1} = []; TempDateTime(cellfun(@(s) isempty(s),TempDateTime))=[];
                
                % If the HarmonicTestFlag is set, the entire analysis is
                % conducted with the frequency of stimulation presumed to
                % be twice its actual value
                
                if (HarmonicTestFlag==1)
                    TempTrialFrequencies=TempTrialFrequencies*2;
                end
                
                % Now remove those trials that have insufficient time points
                % A trial could fail in two ways. It's end point could
                % arrive in time before the minimum trial length specified.
                % Or, it could have an insufficient number of measurements
                % to pass the BadNaNThreshold even before spikes / blinks
                % are accounted for later. Sampling frequency is set to
                % half of the actual frequency given the tendency of
                % initial recordings to capture every-other-measurement.
                
                NumTrials = length(TempData);
                trial = 1;
                while (trial<=NumTrials)
                    if (((TempData(trial).time(end)/1000) < minimum_length_trial)||...
                            (((length(TempData(trial).time))/(sampling_frequency/2)) < (final_trial_length*(1-BadNan))))
                        FirstStageBadTrials=FirstStageBadTrials+1;
                        TempData(trial) = [];
                        TempTrialFrequencies(trial) = [];
                        TempTrialDirections(trial) = [];
                        TempTrialPhases(trial) = [];
                        TempTrialFirstSample(trial) = [];
                        TempTrialFirstSampleInter(trial) = [];
                        TempDateTime{trial} = []; TempDateTime(cellfun(@(s) isempty(s),TempDateTime))=[];
                        NumTrials = NumTrials-1;
                    else trial = trial+1;
                    end
                end
                
                
                % Some experiments have multiple, initial background
                % trials. For now, these are removed. We might
                % figure out something else to do with them later.
                % These background trials are detected by the property
                % of having a trial frequency of zero
                bgAdaptTrial = ~cellfun(@isempty,strfind(TempTrialDirections, 'Background-60s'));
                NumTrials = length(TempData);
                trial = 1;
                while (trial<=NumTrials)
                    if TempTrialFrequencies(trial)==0 && bgAdaptTrial(trial)
                        FirstStageBadTrials=FirstStageBadTrials+1;
                        TempData(trial) = [];
                        TempTrialFrequencies(trial) = [];
                        TempTrialDirections(trial) = [];
                        TempTrialPhases(trial) = [];
                        TempTrialFirstSample(trial) = [];
                        TempTrialFirstSampleInter(trial) = [];
                        TempDateTime{trial} = []; TempDateTime(cellfun(@(s) isempty(s),TempDateTime))=[];
                        NumTrials = NumTrials-1;
                    else trial = trial+1;
                    end
                end
                
                
                
                % Now concatenate TempData to the full Data array
                
                if (FirstGoodSessionFlag==1)
                    Data = TempData;
                    TrialFrequencies = TempTrialFrequencies;
                    TrialDirections = TempTrialDirections;
                    TrialPhases = TempTrialPhases;
                    TrialFirstSampleTime = TempTrialFirstSample;
                    TrialDateTime = TempDateTime;
                    FirstGoodSessionFlag=0;
                else
                    Data(end+1:end+length(TempData)) = TempData(1:end);
                    TrialFrequencies(end+1:end+length(TempTrialFrequencies))=...
                        TempTrialFrequencies(1:end);
                    TrialDirections(end+1:end+length(TempTrialDirections))=...
                        TempTrialDirections(1:end);
                    TrialPhases(end+1:end+length(TempTrialPhases))=...
                        TempTrialPhases(1:end);
                    TrialFirstSampleTime(end+1:end+length(TempTrialFirstSample)) =...
                        TempTrialFirstSample(1:end);
                    TrialDateTime(end+1:end+length(TempDateTime)) = TempDateTime;
                end
            end % if file exists
        end % loop over directions
    end % loop over sessions
    
    fprintf(' Done.'); % Notify user
    fprintf(['\n  - Dropped ' num2str(FirstStageBadTrials) ' trials for insufficient data points or additional background.']); % Update user
    
    % Define some variables to hold the processed data
    
    TimeSeries=zeros(length(Data),datalen);
    MeanPupilByTrial=zeros(1,length(Data));
    AmplitudesByTrial=zeros(1,length(Data));
    PhasesByTrial=zeros(1,length(Data));
    BackgroundAdaptTrials=zeros(length(BackgroundAdaptData),adapt_length*sampling_frequency);
    
    % Smooth and interpolate the adaptation periods and average them
    % across sessions. Also check to see if the ratio interupted measurements
    % for the background period was greater than 0.1 (10%). If so,
    % discard the background trial. Only do this check if the
    % OLFlickerSensitivityFlag is set.
    
    if OLFlickerSensitivityFlag
        
        for session=1:length(BackgroundAdaptData)
            
            if (BackgroundAdaptData(session).ratioInterupt>0.1)
                BackgroundAdaptTrials(session,:)=0;
            else
                
                % Shift the time stamps in the data array back by the
                % StimOnsetDelay value.
                
                BackgroundAdaptData(session).time=BackgroundAdaptData(session).time-StimOnsetDelay;
                
                % smooth, interpolate, and resample the data
                
                iy=SGolaySmooth(BackgroundAdaptData(session).time,BackgroundAdaptData(session).diameter,sgolay_span,sgolay_polynomial,sampling_frequency,adapt_length);
                
                % clip to full_trial_length and store data BackgroundAdaptTrial
                % variable
                
                iy=iy(1:adapt_length*sampling_frequency);
                
                BackgroundAdaptTrials(session,:)=iy;
            end % process background trial given sufficiently small interupt value
        end
    else
        for session=1:length(BackgroundAdaptData)
            
            % Shift the time stamps in the data array back by the
            % StimOnsetDelay value.
            
            BackgroundAdaptData(session).time=BackgroundAdaptData(session).time-StimOnsetDelay;
            
            % smooth, interpolate, and resample the data
            
            iy=SGolaySmooth(BackgroundAdaptData(session).time,BackgroundAdaptData(session).diameter,sgolay_span,sgolay_polynomial,sampling_frequency,adapt_length);
            
            % clip to full_trial_length and store data BackgroundAdaptTrial
            % variable
            
            iy=iy(1:adapt_length*sampling_frequency);
            
            BackgroundAdaptTrials(session,:)=iy;
        end
    end
    
    AverageBackgroundAdapt=nansum(BackgroundAdaptTrials);
    AverageBackgroundAdapt=AverageBackgroundAdapt/length(BackgroundAdaptData);
    
    % Smooth, interpolate, and de-spike the data
    %  Calculate and store mean pupil, amplitude, and phase
    
    fprintf(['\n  - Smoothing and interpolating ' num2str(length(Data)) ' trials:  ']); % Update user
    
    trialCounter=1;
    
    for trial=1:length(Data)
        
        % Run a trial counter unless the TrialInspector flag is set
        
        if (TrialInspectorFlag~=1)
            for j=0:log10(trial-1)
                fprintf('\b');          % delete previous counter display
            end
            fprintf('%d', trial);       % update progress trial number
        else
            fprintf('\n');
        end
        
        % Shift the time stamps in the data array back by the
        % StimOnsetDelay value. This accounts for the delay in the onset of
        % the stimulus modulation relative to the time stamps on the pupil
        % size measurements. This will result in some time stamps having a
        % negative value, but the subsequent SGolaySmooth step will handle
        % this and return a vector that is just from time point zero
        % onward.
        
        Data(trial).time=Data(trial).time-StimOnsetDelay;
        
        % smooth, interpolate, and resample the data
        
        iy=SGolaySmooth(Data(trial).time,Data(trial).diameter,sgolay_span,sgolay_polynomial,sampling_frequency,full_trial_length);
        
        % clip to full_trial_length
        
        iy=iy(1:full_trial_length*sampling_frequency);
        
        % Temporarily crop to final_trial_length to allow corrections of
        % phase offsets in the stimuli.
        % Check to make sure that we are not asking for a final trial
        % length longer than the full trial length
        
        if (final_trial_length>full_trial_length)
            error('final_trial_length must be less than or equal to full_trial_length')
        end
        
        iy=iy((full_trial_length-final_trial_length)*sampling_frequency+1:full_trial_length*sampling_frequency);
        
        
        % Align to zero stimulus phase. Phase shifting of the stimuli is in
        % units of the stimulus modulation, which is typically implemented
        % in 200 discrete steps (Set in AlphaSpacing). The stored trial
        % phase value is the step on which the stimulus started. For a
        % sinusoidal modulation, step 1 corresponds to a a value of zero on
        % a sinusoid at the start of the rising phase. Because phase shifts
        % are in units of stimulus steps, the amount of absolute time shift
        % is proportional to the frequency of stimulation. Aso, the shift
        % is always positive, as the value in TrialPhases is either zero
        % (requiring no shift) or a positive integer.
        
        % One is subtracted from the TrialPhases value as a value of one
        % corresponds to the first mirror setting group, which is a phase
        % advancement of zero.
        
        % This requires that the frequency
        % of stimulation has an intenger number of cycles in the final
        % trial length window. We check for this and issue an error if it
        % is not satisfied
        
        %        if ~isequal((TrialFrequencies(trial)*final_trial_length),floor(TrialFrequencies(trial)*final_trial_length))
        %            error('The stim frequency has a non-intenger number of cycles in the final trial window')
        %        end
        
        % In v2.0 of the pupil experiments, we handle the phase a little
        % differently, namely in terms of actual phase angle (in deg),
        % passed in here as a prop of cycle shifting.
        if OLFlickerSensitivityFlag == 1
            AmountToShift = (TrialPhases(trial))*...  % Proportion of stim cycle phase advanced
                (1/TrialFrequencies(trial))*...  % Length of stimulus cycle in seconds
                sampling_frequency*...;  % number of data samples per second
                1;    % we need to phase advance, therefore -1 here
        else
            AmountToShift= ((TrialPhases(trial)-1)/AlphaSpacing)*...  % Proportion of stim cycle phase advanced
                (1/TrialFrequencies(trial))*...  % Length of stimulus cycle in seconds
                sampling_frequency*...;  % number of data samples per second
                (-1);    % we need to phase advance
        end
        
        % fshift implements a sinc shift, as we may have a non-integer shift
        
        iy=fshift(iy,AmountToShift);
        
        % Restore the time-series to its full length and NaN out the initial
        % values that we do not wish to contribute to the calculation of
        % the amplitude. This allows the
        % phase calculations to proceed on un-shifted data. This is
        % primarily needed for the TTF4D experiment which did not
        % include stimulus phase randomization
        
        if (final_trial_length<full_trial_length)
            iy((full_trial_length-final_trial_length)*sampling_frequency+1:full_trial_length*sampling_frequency)=iy;
            iy(1:(full_trial_length-final_trial_length)*sampling_frequency)=nan;
        end
        
        
        % Mean center the time-series, accounting for the Nan values,
        % and convert to proportion change. Store the mean pupil size.
        
        MeanPupilByTrial(trial)=nansum(iy)/(datalen-sum(isnan(iy)));
        iy=(iy-MeanPupilByTrial(trial))/MeanPupilByTrial(trial);
        
        % Hang on to an original form of iy to compare to the cleaned
        
        orig_iy=iy;
        
        % Include "spike remover" here if we so desire
        
        if (spike_remover_params(1)~=0)
            [iy, indx]=SpikeRemover(iy,spike_remover_params);
        end
        
        % remove those time points that have greater than the BadThreshold
        % percent change.
        removePoints = find(abs(iy)>BadPercentChangeThreshold);
        iy(removePoints)=NaN;
        
        
        
        %%% Analysis of gaze data
        if analyzeGazeData
            
            % resample the gaze data without doing any interpolation
            
            ix = linspace(0,full_trial_length-(1/sampling_frequency),full_trial_length*sampling_frequency);
            iy_gaze_x = interp1(Data(trial).rawTimeStamps/1000,Data(trial).rawMmPositions(:, 1),ix,'spline');  % interpolated data to even sample times, starting from zero
            iy_gaze_y = interp1(Data(trial).rawTimeStamps/1000,Data(trial).rawMmPositions(:, 2),ix,'spline');  % interpolated data to even sample times, starting from zero
            
            
            % We extract the indices of the discarded points from the
            % horizontal position, then remove the corresponding data points as
            % well. These params seem to work.
            [iy_gaze_x, idx] = SpikeRemover(iy_gaze_x,[5 100 3]);
            iy_gaze_y(idx) = NaN;
            
            
            %% Figure out the variance with a sliding window.
            win_width = 120;  %Sliding window width
            slide_incr = 1;  %Slide for each iteration
            numstps = (length(iy_gaze_y)-win_width)/slide_incr; %Number of windows
            mean_win=nan(1,size(iy,2));
            
            for i = 1:numstps
                mean_win(i+floor(win_width/2)) = nanstd(iy_gaze_x(i:i+win_width));  %Calculation for each window
            end
            plot(mean_win)
            
            iy=mean_win;
            
        end
        
        
        % If the TrialInspectorFlag is set, display the time series
        
        if (TrialInspectorFlag==1)
            figure(figTrialInspector);
            plot(orig_iy,'r')
            hold on;
            plot(iy,'k');
            hold off;
            fprintf(['trial: ' int2str(trial) ' ' num2str(TrialFrequencies(trial)) ' ' char(TrialDirections(trial)) ' ' num2str(TrialPhases(trial))]);
        end
        
        % If the time series within the final trial length window is
        % is composed of fewer than BadNanThresh nans then save the trial
        
        if (sum(isnan(iy((full_trial_length-final_trial_length)*sampling_frequency+1:full_trial_length*sampling_frequency)))<...
                (final_trial_length*sampling_frequency*BadNanThreshold))
            
            % If the TrialInspectorFlag is on, tag trial as good
            if (TrialInspectorFlag==1)
                fprintf(' - GOOD');
                pause;
            end
            
            % mean center the time-series again, as we have likely
            % added more Nans
            
            iy=iy-(nansum(iy)/(datalen-sum(isnan(iy))));
            
            % Store the average time series and calculate amplitude and phase
            
            TimeSeries(trialCounter,:)=iy;
            if (StutterErrorFlag==1)
                [AmplitudesByTrial(trialCounter),PhasesByTrial(trialCounter),~] =...
                    LeastSquaresSpectralFit_StutterError(iy,TrialFrequencies(trial),sampling_frequency);
            else % StutterError is false
                [AmplitudesByTrial(trialCounter),PhasesByTrial(trialCounter),~] =...
                    LeastSquaresSpectralFit(iy,TrialFrequencies(trial),sampling_frequency);
            end % StutterErrorFlag if statement
            
            % Copy over the trial values (Freqs, Directions, Phases,
            % MeanPupilSize) to account for any skipped trials due to
            % NaN overage
            
            TrialFrequencies(trialCounter) = TrialFrequencies(trial);
            TrialDirections(trialCounter) = TrialDirections(trial);
            TrialPhases(trialCounter) = TrialPhases(trial);
            MeanPupilByTrial(trialCounter)=MeanPupilByTrial(trial);
            
            % increment the trialCounter
            
            trialCounter=trialCounter+1;
            
        else % check for proportion of NaNs in time series
            if (TrialInspectorFlag==1)
                fprintf(' - BAD');
                pause;
            end
            
        end
        
    end
    
    % Trim the data variables to account for any skipped trials due to
    % exceeding the allowable proportion of NaNs
    
    TimeSeries=TimeSeries(1:trialCounter-1,:);
    MeanPupilByTrial=MeanPupilByTrial(1:trialCounter-1);
    AmplitudesByTrial=AmplitudesByTrial(1:trialCounter-1);
    PhasesByTrial=PhasesByTrial(1:trialCounter-1);
    TrialFrequencies = TrialFrequencies(1:trialCounter-1);
    TrialDirections = TrialDirections(1:trialCounter-1);
    TrialPhases = TrialPhases(1:trialCounter-1);
    TrialDateTime = TrialDateTime(1:trialCounter-1);
    TrialFirstSampleTime = TrialFirstSampleTime(1:trialCounter-1);
    
    if ~(TrialInspectorFlag==1)
        for j=0:log10(trial-1)
            fprintf('\b'); % delete previous counter display
        end
    end
    
    fprintf('Done.'); % notify user we are done the loop
    
    % Save the retention rate in a variable
    nTotalTrials = length(Data);
    nGoodTrials = (trialCounter-1);
    retentionRate = nGoodTrials/nTotalTrials;
    
    % Print out the number of good trials
    fprintf(['\n  - Found ' num2str(trialCounter-1) ' good trials.']);
    fprintf(['\n  - Retention rate is ' num2str(retentionRate*100) '%%.']);
    
    
    % If there are no good trials for the subject, stop at this point
    % so that the operator can remove the subject from analysis
    
    if (trialCounter-1==0)
        fprintf('No good trials found for this subject. Exiting.');
        break
    end
    
    % Allocate variables to hold results
    
    % Create an average time-series for each unique frequency and direction
    % crossing. Calculate for these the best fit, amplitude and phase.
    % Also, calculate a noise level bootstrap for each crossing.
    UniqueFreqs=unique(TrialFrequencies);
    UniqueDirections=unique(TrialDirections);
    
    % A "synthetic" direction is added, which is the sum of LM, S, and  Mel.
    % This is why these results variables are +1 in the Directions dimension
    % when the SyntheticIsoFlag is set to one.
    
    AvgTimeSeries=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag,datalen);
    SEMTimeSeries=zeros(length(UniqueFreqs),length(UniqueDirections),datalen);
    AvgTimeSeriesAmplitudes=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMAmplitudeSplitHalfAvgTimeSeries=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMPhaseSplitHalfAvgTimeSeries=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    AvgTimeSeriesPhases=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    AvgModelFits=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag,datalen);
    NoiseAmplitudes=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMNoiseAmplitudes=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    NoisePhases=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMNoisePhases=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    AvgMeanPupilByFxD=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMMeanPupilByFxD=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    SEMAmplitudesByFxD=zeros(length(UniqueFreqs),length(UniqueDirections)+SyntheticIsoFlag);
    
    % if the HarmonicFlag is set, set up some containers
    % to hold those responses as well
    
    if (HarmonicModelFlag==1)
        AvgTimeSeriesHarmonicAmplitudes=AvgTimeSeriesAmplitudes;
        AvgTimeSeriesHarmonicPhases=AvgTimeSeriesPhases;
    end
    
    fprintf(['\n  - Averaging ' num2str(length(UniqueFreqs)*length(UniqueDirections)) ' crossings of frequency and direction: ']); % Update user
    
    for f=1:length(UniqueFreqs)
        for d=1:length(UniqueDirections)
            IterationCount=(f-1)*length(UniqueDirections) + d;
            for j=0:log10(IterationCount-1)
                fprintf('\b'); % delete previous counter display
            end
            fprintf('%d', IterationCount); % update progress dots
            
            Indices = find(and((TrialFrequencies==UniqueFreqs(f)),strcmp(TrialDirections,UniqueDirections(d))));
            if (not(isempty(Indices)))
                
                % Create the average time series. This is done by first
                %  assembling a matrix across time-series, and then using
                %  the nanmean to handle NaN points
                
                for i=1:length(Indices)
                    if (i==1)
                        TimeSeriesMatrix=[TimeSeries(Indices(i),:)'];
                    else
                        TimeSeriesMatrix=[TimeSeriesMatrix TimeSeries(Indices(i),:)'];
                    end % if the first timeseries
                end % for number of indicies
                
                AvgTimeSeries(f,d,:)=nanmean(TimeSeriesMatrix,2);
                
                TimeSeriesMatrices{f, d} = TimeSeriesMatrix;
                
                % Obtain the SD across the time series. the 'nanstd'
                % command takes a flag that handles normalizing the
                % standard deviation by n or n-1. We go with n
                
                NANSTD_NORMFLAG=1;
                SEMTimeSeries(f,d,:)=nanstd(TimeSeriesMatrix,NANSTD_NORMFLAG,2);
                
                % Now divide by the square root of the number of measures
                % that contributed to each time point in the average time
                % series, taking into account the fact that some time
                % points have nan values in some observations
                
                for t=1:datalen
                    SEMTimeSeries(f,d,t)=SEMTimeSeries(f,d,t)/sqrt(length(Indices)-sum(isnan(TimeSeriesMatrix(t,:))));
                end
                
                % Get the fit to the average temporal respone
                
                if (StutterErrorFlag==1)
                    [AvgTimeSeriesAmplitudes(f,d),AvgTimeSeriesPhases(f,d),AvgModelFits(f,d,:)] =...
                        LeastSquaresSpectralFit_StutterError(transpose(squeeze(AvgTimeSeries(f,d,:))),UniqueFreqs(f),sampling_frequency);
                else
                    [AvgTimeSeriesAmplitudes(f,d),AvgTimeSeriesPhases(f,d),AvgModelFits(f,d,:)] =...
                        LeastSquaresSpectralFit(transpose(squeeze(AvgTimeSeries(f,d,:))),UniqueFreqs(f),sampling_frequency);
                end % StutterErrorFlag if statement
                
                % Check if the HarmonicFlag is set, model the harmonic
                % response at twice the modulation frequency as well
                
                if (HarmonicModelFlag==1)
                    [AvgTimeSeriesHarmonicAmplitudes(f,d),AvgTimeSeriesHarmonicPhases(f,d),HarmonicModelFit] =...
                        LeastSquaresSpectralFit(transpose(squeeze(AvgTimeSeries(f,d,:))),UniqueFreqs(f)*2,sampling_frequency);
                    AvgModelFits(f,d,:)=squeeze(AvgModelFits(f,d,:))+HarmonicModelFit';
                end
                
                AvgMeanPupilByFxD(f,d)=sum(MeanPupilByTrial(Indices))/length(Indices);
                SEMMeanPupilByFxD(f,d)=std(MeanPupilByTrial(Indices))/sqrt(length(Indices));
                
                % The bootstrap analysis returns an estimate of the
                % standard error of the mean
                
                
                [~,~,SEMAmplitudeSplitHalfAvgTimeSeries(f,d),SEMPhaseSplitHalfAvgTimeSeries(f,d)]=...
                    NoiseBootStrap(TimeSeries,Indices,UniqueFreqs(f),sampling_frequency);
                SEMAmplitudeSplitHalfAvgTimeSeries(f,d)=SEMAmplitudeSplitHalfAvgTimeSeries(f,d);
                SEMPhaseSplitHalfAvgTimeSeries(f,d)=SEMPhaseSplitHalfAvgTimeSeries(f,d);
                
                % Get the SEM of amplitude estimates across trials
                
                SEMAmplitudesByFxD(f,d)=std(AmplitudesByTrial(Indices))/sqrt(length(Indices));
                
                % Calculate the mean amplitude of noise at
                % non-stimulated frequencies if there was more than one
                % frequency of stimulation in the experiment.
                
                if (length(UniqueFreqs)>1)
                    NotIndices = find(and((TrialFrequencies~=UniqueFreqs(f)),strcmp(TrialDirections,UniqueDirections(d))));
                    [NoiseAmplitudes(f,d),NoisePhases(f,d),SEMNoiseAmplitudes(f,d),SEMNoisePhases(f,d)]=...
                        NoiseBootStrap(TimeSeries,NotIndices,UniqueFreqs(f),sampling_frequency);
                end
                
            end % if indices is not length zero
        end % for number of unique directions
    end % for number of unique frequencies
    
    for j=0:log10(IterationCount-1)
        fprintf('\b'); % delete previous counter display
    end
    fprintf('Done.'); % notify user we are done the loop
    
    
    % The SyntheticIso is hard-coded to sum the 2nd, 3rd, and 4th directions,
    %  which are LM, Mel, and S. This is ugly and kludgy and will break in
    %  any analysis other than the TTF4D and OPN4
    
    % If the HarmonicTestFlag is also set to one, then calculate the
    % SyntheticIso as the average of the other directions, as we believe
    % that the harmonic is generated after the combination of receptor
    % channels.
    
    % If the OPN4ScaleFlag is set, this indicates that the measured
    % directions included LM, Mel, and S, but did not include Iso itself. A
    % SyntheticIso value is still created, but need to adjust the indexing.
    
    % At times, we may wish to examine other sums / differences of the
    % responses. As a kludge, these sums and differences of modulation
    % directions will still be entered into the "SynIso" array index point.
    % Look at the notes below and flags below!
    
    SumCase = 1; % Create a SynIso from LM+Mel+S
    %    SumCase = 2; % Create Mel+S response and place in the "SynIso" array position
    %    SumCase = 3; % Create LM+S response and place in the "SynIso" array position
    
    if (SyntheticIsoFlag==1)
        if (OPN4ScaleFlag==1)
            for f=1:length(UniqueFreqs)
                if (HarmonicTestFlag==1)
                    if (SumCase==1)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,1,:))+squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:)))/3;
                    end
                    if (SumCase==2)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:)))/2;
                    end
                    if (SumCase==3)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,1,:))+squeeze(AvgTimeSeries(f,3,:)))/2;
                    end
                else
                    if (SumCase==1)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,1,:))+squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:)));
                    end
                    if (SumCase==2)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:)));
                    end
                    if (SumCase==3)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,1,:))+squeeze(AvgTimeSeries(f,3,:)));
                    end
                end
                [AvgTimeSeriesAmplitudes(f,end),AvgTimeSeriesPhases(f,end),AvgModelFits(f,end,:)] = LeastSquaresSpectralFit(transpose(squeeze(AvgTimeSeries(f,end,:))),UniqueFreqs(f),sampling_frequency);
                AvgMeanPupilByFxD(f,end)=(AvgMeanPupilByFxD(f,1)+AvgMeanPupilByFxD(f,2)+AvgMeanPupilByFxD(f,3))/3;
            end % for number of unique frequencies
        else % The OPN4ScaleFlag was not set, so this should be TTF4D data
            for f=1:length(UniqueFreqs)
                if (HarmonicTestFlag==1)
                    if (SumCase==1)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:))+squeeze(AvgTimeSeries(f,4,:)))/3;
                    end
                    if (SumCase==2)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,3,:))+squeeze(AvgTimeSeries(f,4,:)))/2;
                    end
                    if (SumCase==3)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,4,:)))/2;
                    end
                else
                    if (SumCase==1)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,3,:))+squeeze(AvgTimeSeries(f,4,:)));
                    end
                    if (SumCase==2)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,3,:))+squeeze(AvgTimeSeries(f,4,:)));
                    end
                    if (SumCase==3)
                        AvgTimeSeries(f,end,:)=(squeeze(AvgTimeSeries(f,2,:))+squeeze(AvgTimeSeries(f,4,:)));
                    end
                end
                if (StutterErrorFlag==1)
                    [AvgTimeSeriesAmplitudes(f,end),AvgTimeSeriesPhases(f,end),AvgModelFits(f,end,:)] = LeastSquaresSpectralFit_StutterError(transpose(squeeze(AvgTimeSeries(f,end,:))),UniqueFreqs(f),sampling_frequency);
                else
                    [AvgTimeSeriesAmplitudes(f,end),AvgTimeSeriesPhases(f,end),AvgModelFits(f,end,:)] = LeastSquaresSpectralFit(transpose(squeeze(AvgTimeSeries(f,end,:))),UniqueFreqs(f),sampling_frequency);
                end % StutterErrorFlag
                AvgMeanPupilByFxD(f,end)=(AvgMeanPupilByFxD(f,2)+AvgMeanPupilByFxD(f,3)+AvgMeanPupilByFxD(f,4))/3;
            end % for number of unique frequencies
            
        end
    end % if we are creating the synthetic Iso
    
    fprintf('\n  - Plot results.\n'); % Notify user
    
    % Make a different label vector for plotting purposes
    UniqueDirectionLabels = UniqueDirections;
    if (RelabelDirections==1)
        for i = 1:length(oldLabels)
            ix = find(strcmp(UniqueDirectionLabels, oldLabels{i}));
            if ~(isempty(ix))
                UniqueDirectionLabels{ix} = newLabels{i};
            end % a match between the trial label and the replacement label set was found
        end
    end
    if (SyntheticIsoFlag==1)
        if (SumCase==1)
            SynLabel='L/M + Mel + S';
        end
        if (SumCase==2)
            SynLabel='Mel + S';
        end
        if (SumCase==3)
            SynLabel='L/M + S';
        end
        UniqueDirectionLabels=[UniqueDirectionLabels SynLabel];
    end
    
    %% Create the sparkline plots
    
    ErrorMultiplier=2; % Plots ±2SEM on time series sparklines
    
    % first need to determine a constant Y-range across all plot cells
    for f=1:length(UniqueFreqs)
        for d=1:length(UniqueDirections)
            MaxYRange(f,d) = nanmax([ nanmax(squeeze(AvgTimeSeries(f,d,:))+ErrorMultiplier*squeeze(SEMTimeSeries(f,d,:))) nanmax(squeeze(AvgModelFits(f,d,:)))])  ;
            MinYRange(f,d) = nanmin([ nanmin(squeeze(AvgTimeSeries(f,d,:))-ErrorMultiplier*squeeze(SEMTimeSeries(f,d,:))) nanmin(squeeze(AvgModelFits(f,d,:)))])  ;
        end
    end
    yRange=[min(min(MinYRange)) max(max(MaxYRange))];
    yRange(1)=floor(yRange(1)*100)/100;
    yRange(2)=ceil(yRange(2)*100)/100;
    
    % Now make the plots. One column for each direction
    
    figure(figSparkLinePlots);
    SuperPlotTitle=['Subject ' char(SubjectLabel) '/Responses and fits [' num2str(yRange(1)) ',' num2str(yRange(2)) ']'];
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', SuperPlotTitle, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',12)
    if (and(length(UniqueDirections)==1,length(UniqueFreqs)==1))
        PlotTitle=char(UniqueDirectionLabels(1));
        RowLabels=num2str(UniqueFreqs);
        PlotOLPupilSparklines(AvgTimeSeries(:,d,:),AvgModelFits(:,d,:),SEMTimeSeries(:,d,:), ErrorMultiplier, yRange, figSparkLinePlots, RowLabels, PlotTitle, 1, 1);
    else % if we only have one crossing of freq and direction, else...
        
        for d=1:length(UniqueDirections)
            PlotTitle=char(UniqueDirectionLabels(d));
            RowLabels=num2str(UniqueFreqs);
            PlotOLPupilSparklines(squeeze(AvgTimeSeries(:,d,:)),squeeze(AvgModelFits(:,d,:)),squeeze(SEMTimeSeries(:,d,:)),ErrorMultiplier,yRange, figSparkLinePlots, RowLabels, PlotTitle, d, length(UniqueDirections));
        end
    end % if statement catching edge case of one crossing of freq and direction
    
    
    
    %% Create the cycle average sparkline plots
    
    % First need to obtain the cycle average across frequency and
    % direction. Each cycle average is interpolated up to the length of the
    % lowest frequency cycle, to allow all the plots to be displayed on a
    % single graph. Of course, this means that the scaling of the x-axis
    % differs for the plot of each stimulus frequency.
    
    CycleAveragesData=zeros(length(UniqueFreqs),length(UniqueDirections),final_trial_length*sampling_frequency);
    CycleAveragesData(:,:,:)=nan;
    CycleAveragesSEM=zeros(length(UniqueFreqs),length(UniqueDirections),final_trial_length*sampling_frequency);
    CycleAveragesSEM(:,:,:)=nan;
    CycleAveragesModelFit=zeros(length(UniqueFreqs),length(UniqueDirections),final_trial_length*sampling_frequency);
    CycleAveragesModelFit(:,:,:)=nan;
    
    for f=1:length(UniqueFreqs)
        for d=1:length(UniqueDirections)
            
            % Get the cycle average response
            
            MeanCenterFlag=1; % Mean center the data vector before creation of cycle average
            TmpCycleAverage=MakeCycleAverage(squeeze(squeeze(AvgTimeSeries(f,d,:))),sampling_frequency,UniqueFreqs(f),MeanCenterFlag);
            TmpCycleAverage=TmpCycleAverage-nanmean(TmpCycleAverage);
            %CycleAveragesData(f,d,1:length(TmpCycleAverage))=TmpCycleAverage;
            
            % Handle the kludge here that we return two cycles if the cycle
            % length is less thatn 200 samples
            
            FitFrequency=1/final_trial_length;
            if (length(TmpCycleAverage)<200)
                FitFrequency=FitFrequency*2;
            end
            
            % Interpolate up the response
            
            if (length(TmpCycleAverage)==(final_trial_length*sampling_frequency))
                CycleAveragesData(f,d,:)=TmpCycleAverage;
            else
                CycleAveragesData(f,d,:)=sinc_interp(TmpCycleAverage',linspace(0,1,length(TmpCycleAverage)),linspace(0,1,final_trial_length*sampling_frequency));
            end
            
            % Fit a sinusoid to the interpolated cycle average data
            
            [tmpa,tmpb,TmpBestFit] = LeastSquaresSpectralFit(squeeze(squeeze(CycleAveragesData(f,d,:)))',FitFrequency,sampling_frequency);
            
            % Check if the HarmonicFlag is set, model the harmonic
            % response at twice the modulation frequency as well
            
            if (HarmonicModelFlag==1)
                [tmpa,tmpb,TmpBestFitHarmonic] = LeastSquaresSpectralFit(squeeze(squeeze(CycleAveragesData(f,d,:)))',FitFrequency*2,sampling_frequency);
                TmpBestFit=TmpBestFit+TmpBestFitHarmonic;
            end
            
            CycleAveragesModelFit(f,d,:)=TmpBestFit;
            
            % Grab the standard error of the mean of the trial time points
            
            MeanCenterFlag=0; % Do not mean center the SEM vector
            TmpCycleAverage=MakeCycleAverage(squeeze(squeeze(SEMTimeSeries(f,d,:))),sampling_frequency,UniqueFreqs(f),MeanCenterFlag);
            if (length(TmpCycleAverage)==(final_trial_length*sampling_frequency))
                CycleAveragesSEM(f,d,:)=TmpCycleAverage;
            else
                CycleAveragesSEM(f,d,:)=sinc_interp(TmpCycleAverage',linspace(0,1,length(TmpCycleAverage)),linspace(0,1,final_trial_length*sampling_frequency));
            end
            
        end
    end
    
    % first need to determine a constant Y-range across all plot cells
    for f=1:length(UniqueFreqs)
        for d=1:length(UniqueDirections)
            MaxYRange(f,d) = nanmax([ nanmax(squeeze(CycleAveragesData(f,d,:))+ErrorMultiplier*squeeze(CycleAveragesSEM(f,d,:))) nanmax(squeeze(CycleAveragesModelFit(f,d,:)))])  ;
            MinYRange(f,d) = nanmin([ nanmin(squeeze(CycleAveragesData(f,d,:))-ErrorMultiplier*squeeze(CycleAveragesSEM(f,d,:))) nanmin(squeeze(CycleAveragesModelFit(f,d,:)))])  ;
        end
    end
    yRange=[min(min(MinYRange)) max(max(MaxYRange))];
    yRange(1)=floor(yRange(1)*100)/100;
    yRange(2)=ceil(yRange(2)*100)/100;
    
    % Now make the plots. One column for each direction
    
    figure(figCycleAveragePlots);
    SuperPlotTitle=['Subject ' char(SubjectLabel) '/Responses and fits [' num2str(yRange(1)) ',' num2str(yRange(2)) ']'];
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', SuperPlotTitle, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',12)
    if (and(length(UniqueDirections)==1,length(UniqueFreqs)==1))
        PlotTitle=char(UniqueDirectionLabels(1));
        RowLabels=num2str(UniqueFreqs);
        PlotOLPupilSparklines(CycleAveragesData(:,d,:),CycleAveragesModelFit(:,d,:),CycleAveragesSEM(:,d,:), ErrorMultiplier, yRange, figCycleAveragePlots, RowLabels, PlotTitle, 1, 1);
    else % if we only have one crossing of freq and direction, else...
        
        for d=1:length(UniqueDirections)
            PlotTitle=char(UniqueDirectionLabels(d));
            RowLabels=num2str(UniqueFreqs);
            PlotOLPupilSparklines(squeeze(CycleAveragesData(:,d,:)),squeeze(CycleAveragesModelFit(:,d,:)),squeeze(CycleAveragesSEM(:,d,:)),ErrorMultiplier,yRange, figCycleAveragePlots, RowLabels, PlotTitle, d, length(UniqueDirections));
        end
    end % if statement catching edge case of one crossing of freq and direction
    
    %% Determine if each subject is in a separate subplot, or if all the
    %% subject results are stacked
    
    SubjectIDPass=SubjectID;
    NumSubjectsPass=length(Subjects);
    if (StackPlotFlag==1)
        SubjectIDPass=1;
        NumSubjectsPass=1;
    end
    
    %% Plot background adaptation traces
    
    PlotData=AverageBackgroundAdapt(1:end);
    yRange=[0 10];
    figure(figBackgroundAdapt);
    SuperPlotTitle=['10s adaptation to background [0 10mm]'];
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', SuperPlotTitle, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',12)
    PlotTitle=[char(SubjectLabel)];
    RowLabels=['adapt'];
    PlotOLPupilSparklines(PlotData,PlotData*0,PlotData*0,1, yRange,figBackgroundAdapt, RowLabels, PlotTitle, SubjectIDPass, NumSubjectsPass);
    
    
    %% Plot avg pupil size by direction and frequency
    
    PlotData=AvgMeanPupilByFxD;
    PlotErrors=SEMMeanPupilByFxD*2;
    
    PlotTitle=['Subject ' char(SubjectLabel) '/Mean Pupil Diameter [±2SEM]'];
    PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,PlotData,PlotErrors,PlotTitle,'Pupil diameter [mm]',[0 max(max(PlotData))+max(max(PlotErrors))],figMeanPupilByDirection,SubjectIDPass,NumSubjectsPass);
    
    
    %% Plot channel noise by direction
    
    PlotData=NoiseAmplitudes*100;
    PlotErrors=SEMNoiseAmplitudes*2*100;
    
    PlotTitle=['Subject ' char(SubjectLabel) '/Channel noise by direction [±2SEM]'];
    PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,PlotData,PlotErrors,PlotTitle,'Amplitude [% change]',[-0.5 10],figNoiseByDirection,SubjectIDPass,NumSubjectsPass);
    
    
    %% Plot phase by direction
    PlotData=AvgTimeSeriesPhases;
    PlotErrors=SEMPhaseSplitHalfAvgTimeSeries;
    
    PlotTitle=['Subject ' char(SubjectLabel) '/Phase [±2SEM]'];
    PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,PlotData,PlotErrors*2,PlotTitle,'Phase [rads]',[-pi pi],figPhase,SubjectIDPass,NumSubjectsPass);
    
    
    %% Add the channel noise measurement
    UniqueDirectionLabels=[UniqueDirectionLabels 'Noise'];
    AvgTimeSeriesAmplitudes = [AvgTimeSeriesAmplitudes mean(NoiseAmplitudes(:,1:end-1),2)];
    SEMAmplitudeSplitHalfAvgTimeSeries=[SEMAmplitudeSplitHalfAvgTimeSeries mean(SEMNoiseAmplitudes(:,1:end-1),2)];
    AvgTimeSeriesPhases = [AvgTimeSeriesPhases mean(NoisePhases(:,1:end-1),2)];
    SEMPhaseSplitHalfAvgTimeSeries=[SEMPhaseSplitHalfAvgTimeSeries mean(SEMNoisePhases(:,1:end-1),2)];
    
    
    %% Plot amplitude by direction
    PlotData=AvgTimeSeriesAmplitudes*100;
    PlotErrors=SEMAmplitudeSplitHalfAvgTimeSeries*2*100;
    
    PlotTitle=['Subject ' char(SubjectLabel) '/Amplitude [±2SEM]'];
    PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,PlotData,PlotErrors,PlotTitle,'Amplitude [% change]',[-0.5 10],figAmplitude,SubjectIDPass,NumSubjectsPass);
    
    %% Create polar plots
    for f=1:length(UniqueFreqs)
        PlotTitle=['Subject ' char(SubjectLabel) '/' num2str(UniqueFreqs(f)) ' Hz'];
        PlotOLPupilPolar(AvgTimeSeriesAmplitudes(f,:)*100,AvgTimeSeriesPhases(f,:),SEMAmplitudeSplitHalfAvgTimeSeries(f,:)*2*100,SEMPhaseSplitHalfAvgTimeSeries(f,:)*2,max(max(AvgTimeSeriesAmplitudes))*100,PlotTitle,figPolarPlots,f,length(UniqueFreqs),UniqueDirectionLabels);
    end
    
    %% If there is a SyntheticIsoFlag set, create a set of polar plots in
    % which all SynIso responses are scaled to have an amplitude of
    % 1 and a phase of zero.
    
    % Note that the SumCase flag (which is set earlier in this routine at
    % the time of the calculation of SynIso) controls what actually ends up
    % in the "SynIso" array position.
    
    if (SyntheticIsoFlag==1)
        AlignedPhases=AvgTimeSeriesPhases;
        ScaledAmplitudes=AvgTimeSeriesAmplitudes;
        ScaledSEMAmplitudeSplitHalfAvgTimeSeries=SEMAmplitudeSplitHalfAvgTimeSeries;
        if (SumCase==1)
            MaxAmp=1.25;
        end
        if (SumCase==2)
            MaxAmp=30;
        end
        if (SumCase==3)
            MaxAmp=3;
        end
        
        if (OPN4ScaleFlag==1)
            SynIsoIndex=4;
            for f=1:length(UniqueFreqs)
                AlignedPhases(f,:)=AvgTimeSeriesPhases(f,:)-AvgTimeSeriesPhases(f,SynIsoIndex);
                IsoAmp=AvgTimeSeriesAmplitudes(f,SynIsoIndex);
                ScaledAmplitudes(f,:)=AvgTimeSeriesAmplitudes(f,:)/IsoAmp;
                ScaledSEMAmplitudeSplitHalfAvgTimeSeries(f,:)=SEMAmplitudeSplitHalfAvgTimeSeries(f,:)/IsoAmp;
                PlotTitle=['Subject ' char(SubjectLabel) '/' num2str(UniqueFreqs(f)) ' Hz'];
                PlotOLPupilPolar(ScaledAmplitudes(f,:),AlignedPhases(f,:),ScaledSEMAmplitudeSplitHalfAvgTimeSeries(f,:)*2,SEMPhaseSplitHalfAvgTimeSeries(f,:)*2,MaxAmp,PlotTitle,figScaledPolarPlots,f,length(UniqueFreqs),UniqueDirectionLabels);
            end
        else
            SynIsoIndex=5;
            for f=1:length(UniqueFreqs)
                AlignedPhases(f,:)=AvgTimeSeriesPhases(f,:)-AvgTimeSeriesPhases(f,SynIsoIndex);
                IsoAmp=AvgTimeSeriesAmplitudes(f,SynIsoIndex);
                ScaledAmplitudes(f,:)=AvgTimeSeriesAmplitudes(f,:)/IsoAmp;
                ScaledSEMAmplitudeSplitHalfAvgTimeSeries(f,:)=SEMAmplitudeSplitHalfAvgTimeSeries(f,:)/IsoAmp;
                PlotTitle=['Subject ' char(SubjectLabel) '/' num2str(UniqueFreqs(f)) ' Hz'];
                PlotOLPupilPolar(ScaledAmplitudes(f,:),AlignedPhases(f,:),ScaledSEMAmplitudeSplitHalfAvgTimeSeries(f,:)*2,SEMPhaseSplitHalfAvgTimeSeries(f,:)*2,MaxAmp,PlotTitle,figScaledPolarPlots,f,length(UniqueFreqs),UniqueDirectionLabels);
            end
        end
    end
    
    
    %% Save out the data
    
    if (SaveDataFlag==1)
        
        currDir=pwd;
        
        if ~isdir(ResultsFullPath)
            mkdir(ResultsFullPath)
        end
        cd(ResultsFullPath)
        
        cd(ResultsFullPath)
        if ~isdir(ResultsFullPath)
            mkdir(ResultsFullPath)
        end
        
        if ~isdir(fullfile(ResultsFullPath, ResultsDirName))
            mkdir(ResultsDirName)
        end
        
        % Make a folder for the diagnostic output
        cd(ResultsDirName);
        diagnosticDirName = 'Diagnostic';
        if ~isdir(fullfile(ResultsFullPath, diagnosticDirName))
            mkdir(diagnosticDirName)
        end
        
        dataFile = sprintf('%s-results.csv', Subjects{SubjectID});
        c = CSVFile(dataFile, true);
        
        phase = AvgTimeSeriesPhases;
        amplitude = AvgTimeSeriesAmplitudes;
        phaseError = SEMPhaseSplitHalfAvgTimeSeries;
        amplitudeError = SEMAmplitudeSplitHalfAvgTimeSeries;
        
        allDirections = repmat(UniqueDirectionLabels, length(UniqueFreqs), 1);
        c = c.addColumn('Direction', 's');
        c = c.setColumnData('Direction', allDirections(:));
        
        c = c.addColumn('Frequency [Hz]', 'g');
        c = c.setColumnData('Frequency [Hz]', repmat(UniqueFreqs', length(UniqueDirectionLabels), 1));
        
        c = c.addColumn('Pupil size mean [mm]', 'g');
        c = c.setColumnData('Pupil size mean [mm]', AvgMeanPupilByFxD(:));
        
        c = c.addColumn('Pupil size error [1SEM]', 'g');
        c = c.setColumnData('Pupil size error [1SEM]', SEMMeanPupilByFxD(:));
        
        c = c.addColumn('Amplitude mean [proportion change]', 'g');
        c = c.setColumnData('Amplitude mean [proportion change]', amplitude(:));
        
        c = c.addColumn('Phase mean (±pi)', 'g');
        c = c.setColumnData('Phase mean (±pi)', phase(:));
        
        c = c.addColumn('Amplitude error [1SEM]', 'g');
        c = c.setColumnData('Amplitude error [1SEM]', amplitudeError(:));
        
        c = c.addColumn('Phase error [1SEM]', 'g');
        c = c.setColumnData('Phase error [1SEM]', phaseError(:));
        
        % Add contrasts (if config filenames were passed)
        if ~(isempty(configFileNames))
            for r = 1:length(receptorNames)
                % Target contrasts
                c = c.addColumn([receptorNames{r} ' contrast (target)'], 'g');
                tmpData = (round(repmat(contrastTarget(:, r), 1, length(UniqueFreqs))*1000)/1000)';
                c = c.setColumnData([receptorNames{r} ' contrast (target)'], tmpData(:));
                
                % Predicted contasts (with dark)
                c = c.addColumn([receptorNames{r} ' contrast (predicted)'], 'g');
                tmpData = (round(repmat(contrastPredicted(:, r), 1, length(UniqueFreqs))*1000)/1000)';
                c = c.setColumnData([receptorNames{r} ' contrast (predicted)'], tmpData(:));
            end
        end
        
        c.write;
        
        % Save out information about data quality: # trials, # good trials,
        % retention rate
        diagnosticFileName = sprintf('%s-data_quality.csv', Subjects{SubjectID});;
        fid = fopen(fullfile(ResultsFullPath, ResultsDirName, diagnosticDirName, diagnosticFileName), 'w');
        fprintf(fid, '%s,%s,%s', '# trials', '# good trials', 'Retention rate');
        fprintf(fid, '\n');
        fprintf(fid, '%g,%g,%g', nTotalTrials, nGoodTrials, retentionRate);
        fclose(fid);
        
        %% Save out the cycle averages
        
        signalAvgCycle.theFreqs = UniqueFreqs;
        signalAvgCycle.theDirections = UniqueDirections;
        signalAvgCycle.theData = CycleAveragesData;
        signalAvgCycle.theDataSEM = CycleAveragesSEM;
        signalAvgCycle.theFit = CycleAveragesModelFit;
        
        if HarmonicModelFlag
            save([Subjects{SubjectID} '-cycleaverageharmonic'], 'signalAvgCycle');
        else
            save([Subjects{SubjectID} '-cycleaverage'], 'signalAvgCycle');
        end
        
        cd(currDir);
    end % if SaveDataFlag
    
    
    
    %% If there is a SyntheticIsoFlag set, write out the data in
    % which all SynIso responses are scaled to have an amplitude of
    % 1 and a phase of zero.
    
    if and((SaveDataFlag==1),(SyntheticIsoFlag==1))
        
        currDir=pwd;
        
        
        if ~isdir(ResultsFullPath)
            mkdir(ResultsFullPath)
        end
        cd(ResultsFullPath)
        
        if ~isdir(fullfile(ResultsFullPath, ResultsDirName))
            mkdir(ResultsDirName)
        end
        
        % Make a folder for the diagnostic output
        cd(ResultsDirName);
        diagnosticDirName = 'Diagnostic';
        if ~isdir(fullfile(ResultsFullPath, diagnosticDirName))
            mkdir(diagnosticDirName)
        end
        
        if (SumCase==1)
            dataFile = sprintf('%s-IsoScaled_results.csv', Subjects{SubjectID});
        end
        if (SumCase==2)
            dataFile = sprintf('%s-Mel+S_Scaled_results.csv', Subjects{SubjectID});
        end
        if (SumCase==3)
            dataFile = sprintf('%s-L+M+S_Scaled_results.csv', Subjects{SubjectID});
        end
        
        c = CSVFile(dataFile, true);
        
        phase = AlignedPhases;
        amplitude = ScaledAmplitudes;
        phaseError = SEMPhaseSplitHalfAvgTimeSeries;
        amplitudeError = ScaledSEMAmplitudeSplitHalfAvgTimeSeries;
        
        allDirections = repmat(UniqueDirectionLabels, length(UniqueFreqs), 1);
        c = c.addColumn('Direction', 's');
        c = c.setColumnData('Direction', allDirections(:));
        
        c = c.addColumn('Frequency [Hz]', 'g');
        c = c.setColumnData('Frequency [Hz]', repmat(UniqueFreqs', length(UniqueDirectionLabels), 1));
        
        c = c.addColumn('Pupil size mean [mm]', 'g');
        c = c.setColumnData('Pupil size mean [mm]', AvgMeanPupilByFxD(:));
        
        c = c.addColumn('Pupil size error [1SEM]', 'g');
        c = c.setColumnData('Pupil size error [1SEM]', SEMMeanPupilByFxD(:));
        
        c = c.addColumn('Amplitude mean [proportion change]', 'g');
        c = c.setColumnData('Amplitude mean [proportion change]', amplitude(:));
        
        c = c.addColumn('Phase mean (±pi)', 'g');
        c = c.setColumnData('Phase mean (±pi)', phase(:));
        
        c = c.addColumn('Amplitude error [1SEM]', 'g');
        c = c.setColumnData('Amplitude error [1SEM]', amplitudeError(:));
        
        c = c.addColumn('Phase error [1SEM]', 'g');
        c = c.setColumnData('Phase error [1SEM]', phaseError(:));
        
        % Add contrasts (if config filenames were passed)
        if ~(isempty(configFileNames))
            for r = 1:length(receptorNames)
                % Target contrasts
                c = c.addColumn([receptorNames{r} ' contrast (target)'], 'g');
                tmpData = (round(repmat(contrastTarget(:, r), 1, length(UniqueFreqs))*1000)/1000)';
                c = c.setColumnData([receptorNames{r} ' contrast (target)'], tmpData(:));
                
                % Predicted contasts (with dark)
                c = c.addColumn([receptorNames{r} ' contrast (predicted)'], 'g');
                tmpData = (round(repmat(contrastPredicted(:, r), 1, length(UniqueFreqs))*1000)/1000)';
                c = c.setColumnData([receptorNames{r} ' contrast (predicted)'], tmpData(:));
            end
        end
        
        c.write;
        
        % Save out information about data quality: # trials, # good trials,
        % retention rate
        diagnosticFileName = sprintf('%s-data_quality.csv', Subjects{SubjectID});;
        fid = fopen(fullfile(ResultsFullPath, ResultsDirName, diagnosticDirName, diagnosticFileName), 'w');
        fprintf(fid, '%s,%s,%s', '# trials', '# good trials', 'Retention rate');
        fprintf(fid, '\n');
        fprintf(fid, '%g,%g,%g', nTotalTrials, nGoodTrials, retentionRate);
        fclose(fid);
        
        cd(currDir);
    end % if SaveDataFlag and SyntheticIso
    
    if (~isempty(cell2mat(strfind(Protocols, 'Luxotonic')))) || ~isempty(cell2mat(strfind(Protocols, 'Wave')))
        AvgTimeSeries = squeeze(AvgTimeSeries);
        
        AvgIso = AvgTimeSeries(1, :);
        AvgMel = AvgTimeSeries(2, :);
        AvgIso(1:480) = [];
        AvgMel(1:480) = [];
        AvgIsoCyc = reshape(AvgIso, 960, 6);
        AvgMelCyc = reshape(AvgMel, 960, 6);
            

        luxFig = figure;
        subplot(1, 2, 1);
        t = linspace(0, 48, 960);
        a = shadedErrorBar(t, mean(AvgIsoCyc, 2), 1*std(AvgIsoCyc, [], 2)/sqrt(6), {'LineWidth', 2, 'Color', 'k'}); hold on;
        %b = plot(t, mean(AvgLMSCyc, 2)+mean(AvgMelCyc, 2), '--r', 'LineWidth', 2);
        plot(t, (7.5/100)+0.009*square(2*pi*1/48*t-pi), '-k')
        plot([0 48], [0 0], '--k');
        pbaspect([1 1 1]); ylim([-0.15 0.15]);
        xlabel('t [s]'); ylabel('Pupil change [%]');
        set(gca, 'XTick', [0 24 48]); set(gca, 'YTick', [-0.08 -0.04 0 0.04 0.08], 'YTickLabel', [-8 -4 0 4 8]);
        xlim([-1 49]);
        legend([ a.patch],  '\pm1SEM', 'Location', 'SouthWest'); legend boxoff;  box off;
        title('LMS');
        % Inset
       
        subplot(1, 2, 2);
        a = shadedErrorBar(t, mean(AvgMelCyc, 2), 1*std(AvgMelCyc, [], 2)/sqrt(6), {'LineWidth', 2, 'Color', 'k'}); hold on;
        plot(t, (7.5/100)+0.009*square(2*pi*1/48*t-pi), '-k')
        plot([0 48], [0 0], '--k');
        pbaspect([1 1 1]); ylim([-0.15 0.15]);
        xlabel('t [s]'); ylabel('Pupil change [%]');
        set(gca, 'XTick', [0 24 48]); set(gca, 'YTick', [-0.08 -0.04 0 0.04 0.08], 'YTickLabel', [-8 -4 0 4 8]);
        xlim([-1 49]);
        legend([ a.patch], '\pm1SEM', 'Location', 'SouthWest'); legend boxoff; box off;
        title('Melanopsin');
        set(luxFig, 'PaperPosition', [0 0 20 12])
        set(luxFig, 'PaperSize', [20 12]); %Set the paper to have width 5 and height 5.
        saveas(luxFig, [char(SubjectLabel) '_cycleAvg'], 'pdf');
    end
    
    % Save subject-specific plots
    
    if (SavePlotsFlag==1)
        
        fprintf('  - Now saving plots to disk.\n'); % Notify user
        
        currDir=pwd;
        
        
        if ~isdir(ResultsFullPath)
            mkdir(ResultsFullPath)
        end
        cd(ResultsFullPath)
        
        if ~isdir(fullfile(ResultsFullPath, ResultsDirName))
            mkdir(ResultsDirName)
        end
        cd(ResultsDirName);
        
        % Put everything in the 'Plots' folder
        plotDirName = 'Plots';
        if ~isdir(fullfile(ResultsFullPath, ResultsDirName, plotDirName))
            mkdir(plotDirName)
        end
        cd(plotDirName)
        
        if (RelabelSubjectsFlag==1)
            SubjectLabel = ['S' num2str(SubjectID)];
        else
            SubjectLabel=Subjects(SubjectID);
        end
        
        set(figPolarPlots, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figPolarPlots, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figPolarPlots, [char(SubjectLabel) '-polar_plot.pdf'], 'pdf');
        
        set(figSparkLinePlots, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figSparkLinePlots, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figSparkLinePlots, [char(SubjectLabel) '-sparkline_plot.pdf'], 'pdf');
        
        set(figCycleAveragePlots, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figCycleAveragePlots, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figCycleAveragePlots, [char(SubjectLabel) '-cycleaverage_plot.pdf'], 'pdf');
        
        if (SyntheticIsoFlag==1)
            set(figScaledPolarPlots, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
            set(figScaledPolarPlots, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
            
            if (SumCase==1)
                saveas(figScaledPolarPlots, [char(SubjectLabel) '-IsoScaled_polar_plot.pdf'], 'pdf');
            end
            if (SumCase==2)
                saveas(figScaledPolarPlots, [char(SubjectLabel) '-Mel+S_Scaled_polar_plot.pdf'], 'pdf');
            end
            if (SumCase==3)
                saveas(figScaledPolarPlots, [char(SubjectLabel) '-L+M+S_Scaled_polar_plot.pdf'], 'pdf');
            end
            
        end
        
        % Save channel noise figure
        set(figNoiseByDirection, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figNoiseByDirection, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figNoiseByDirection, [char(SubjectLabel) '-channel_noise.pdf'], 'pdf');
        
        % Save Avg Pupil size figure
        set(figMeanPupilByDirection, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figMeanPupilByDirection, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figMeanPupilByDirection, [char(SubjectLabel) '-mean_pupil_size.pdf'], 'pdf');
        
        % Save amplitude figure
        set(figAmplitude, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figAmplitude, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figAmplitude, [char(SubjectLabel) '-amplitude.pdf'], 'pdf');
        
        % Save phase figure
        set(figPhase, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figPhase, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figPhase, [char(SubjectLabel) '-phase.pdf'], 'pdf');
        
        % Save background adaptation figure
        set(figBackgroundAdapt, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
        set(figBackgroundAdapt, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
        saveas(figBackgroundAdapt, [char(SubjectLabel) '-background_adapt.pdf'], 'pdf');
        
        
        % This function creates the average power spectrum of the TTF4D
        % data after removing effects at odd and even harmonics of the
        % stimulation frequency. The line is commented out and only
        % should be enabled if one wishes to re-run this control
        % analysis that was needed for the PNAS paper.
        
        %[AvgPS]=MakeAvgPS(TimeSeries,TrialFrequencies,TrialDirections,UniqueDirectionLabels,sampling_frequency,SubjectLabel);
        
        
        
        cd(currDir);
        
    end  % if SavePlotsFlag
    
    
    % Wrap stauff around
    
end % subject loop

fprintf('\n*** Done all subjects.\n'); % Notify user

% break here

keyboard;

end % main
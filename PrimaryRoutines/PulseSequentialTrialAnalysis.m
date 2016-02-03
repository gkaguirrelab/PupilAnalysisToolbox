function PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath)
% PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ROUTINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop across subjects and perform the analysis
for SubjectID=1:length(Subjects)
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
    clear MeanMatrix;
    clear AvgTimeSeries;
    clear SEMTimeSeries;
    clear TimeSeriesMatrixStore;
    clear Indices;
    clear TimeSeriesMatrix;
    clear theValid;
    clear theMeanMatrix;
    
    % Prepare plot figures for this subject
    if (params.TrialInspectorFlag==1)
        figTrialInspector=figure;
        pause on;
    end
    
    % Display some information
    fprintf(['\n>>> Processing subject '...
        char(Subjects(SubjectID))]); % Notify user
    fprintf(['\n\t*** Loading and concatenating sessions for subject '...
        char(Subjects(SubjectID))]); % Notify user
    
    % Loop over sessions and directions for a subject and concatenate the data
    FirstGoodSessionFlag=1;
    for sess=1:3
        for p=1:length(Protocols)
            % Load data set
            [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempAlphaSpacing,TempDateTime]=...
                PupilAnalysisToolbox_LoadData(basePath, char(Protocols(p)),char(Subjects(SubjectID)),sess);
            
            if (not(isempty(TempData)))
                fprintf('.'); % Notify user
                
                % The first two trials are discarded from the primary analysis for each session, as they
                % correspond to the adaptation period and the initial wrap-around trial
                % The first trial is kept, however, and averaged for the
                % subject across all directions to provide the response for
                % adaptation to the background.
                BackgroundAdaptData(sess,p)=TempData(1);
                
                TempData(1)=[]; TempTrialFrequencies(1)=[]; TempTrialDirections(1)=[]; TempTrialPhases(1)=[];
                TempData(2)=[]; TempTrialFrequencies(2)=[]; TempTrialDirections(2)=[]; TempTrialPhases(2)=[];
                
                % Now remove those trials that have insufficient time points
                NumTrials = length(TempData);
                trial = 1;
                while (trial<NumTrials)
                    if ((TempData(trial).time(end)/1000) < params.minimum_length_trial)
                        TempData(trial) = [];
                        TempTrialFrequencies(trial) = [];
                        TempTrialDirections(trial) = [];
                        TempTrialPhases(trial) = [];
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
                    FirstGoodSessionFlag=0;
                else
                    Data(end+1:end+length(TempData)) = TempData(1:end);
                    TrialFrequencies(end+1:end+length(TempTrialFrequencies))=...
                        TempTrialFrequencies(1:end);
                    TrialDirections(end+1:end+length(TempTrialDirections))=...
                        TempTrialDirections(1:end);
                    TrialPhases(end+1:end+length(TempTrialPhases))=...
                        TempTrialPhases(1:end);
                end
            end % if file exists
        end % loop over directions
    end % loop over sessions
    
    fprintf('done.'); % Notify user
    
    % Smooth, interpolate, and de-spike the data
    fprintf(['\n\t*** Smoothing and interpolating ' num2str(length(Data)) ' trials:  ']); % Update user
    trialCounter=1;
    
    for trial=1:length(Data)
        % Run a trial counter unless the TrialInspector flag is set
        if (params.TrialInspectorFlag~=1)
            for j=0:log10(trial-1)
                fprintf('\b');          % delete previous counter display
            end
            fprintf('%d', trial);       % update progress trial number
        else
            fprintf('\n');
        end
        
        % Shift the time stamps in the data array back by the
        % params.StimOnsetDelay value. This accounts for the delay in the onset of
        % the stimulus modulation relative to the time stamps on the pupil
        % size measurements. This will result in some time stamps having a
        % negative value, but the subsequent SGolaySmooth step will handle
        % this and return a vector that is just from time point zero
        % onward.
        Data(trial).time=Data(trial).time-params.StimOnsetDelay;
        
        % smooth, interpolate, and resample the data
        iy=SGolaySmooth(Data(trial).time,Data(trial).diameter,params.sgolay_span,params.sgolay_polynomial,params.sampling_frequency,params.full_trial_length);
        
        % clip to params.full_trial_length
        iy=iy(1:params.full_trial_length*params.sampling_frequency);
        
        % Temporarily crop to params.final_trial_length to allow corrections of
        % phase offsets in the stimuli.
        % Check to make sure that we are not asking for a final trial
        % length longer than the full trial length
        if (params.final_trial_length>params.full_trial_length)
            error('params.final_trial_length must be less than or equal to params.full_trial_length')
        end
        iy=iy((params.full_trial_length-params.final_trial_length)*params.sampling_frequency+1:params.full_trial_length*params.sampling_frequency);
        
        % Align to zero stimulus phase. Phase shifting of the stimuli is in
        % units of the stimulus modulation, which is typically implemented
        % in 200 discrete steps (Set in AlphaSpacing). The stored trial
        % phase value is the step on which the stimulus started. For the
        % Gaussian pulse stimuli, the phase shift is performed simply as a
        % shift in the specified mirror settings. The shift direction
        % is always positive, as the value in TrialPhases is either zero
        % (requiring no shift) or a positive integer.
        
        % One is subtracted from the TrialPhases value as a value of one
        % corresponds to the first mirror setting group, which is a phase
        % advancement of zero.
        
        AmountToShift= TrialPhases(trial)*...  % Length of stimulus cycle in seconds
            params.sampling_frequency*...;  % number of data samples per second
            (1);    % we need to phase advance
        
        % fshift implements a sinc shift, as we may have a non-integer shift
        orig_iy = iy;
        %iy=fshift(orig_iy,AmountToShift);
        tmp = NaN*ones(size(orig_iy));
        tmp(1:length((AmountToShift+1):length(iy))) = orig_iy((AmountToShift+1):length(iy));
        TimeSeries(trial, :) = tmp;
        % Store the mean pupil size.
        
        if (params.TrialInspectorFlag==1)
            figure(figTrialInspector);
            plot(iy,'k');
            hold off;
            fprintf(['trial: ' int2str(trial) ' ' num2str(TrialFrequencies(trial)) ' ' char(TrialDirections(trial)) ' ' num2str(TrialPhases(trial))]);
            pause;
        end
        
    end
    
    trialCounter = trial+1;
    
    % Trim the data variables to account for any skipped trials due to
    % exceeding the allowable proportion of NaNs
    TimeSeries=TimeSeries(1:trialCounter-1,:);
    TrialFrequencies = TrialFrequencies(1:trialCounter-1);
    TrialDirections = TrialDirections(1:trialCounter-1);
    TrialPhases = TrialPhases(1:trialCounter-1);
    
    if ~(params.TrialInspectorFlag==1)
        for j=0:log10(trial-1)
            fprintf('\b'); % delete previous counter display
        end
    end
    
    fprintf('done.'); % notify user we are done the loop
    fprintf(['\n\t*** Found ' num2str(trialCounter-1) ' good trials.']);
    
    % Allocate variables to hold results
    % Create an average time-series for each unique frequency and direction
    % crossing. Calculate for these the best fit, amplitude and phase.
    % Also, calculate a noise level bootstrap for each crossing.
    UniqueFreqs=unique(TrialFrequencies);
    UniqueDirections=unique(TrialDirections);
    
    % Clean up the direction labels
    for d = 1:length(UniqueDirections)
        tmp0 = allwords(UniqueDirections{d}, '-');
        UniqueDirectionLabels{d} = [tmp0{2} '-' tmp0{3}];
    end
    
    if (params.RelabelDirections==1)
        for i = 1:length(oldLabels)
            ix = find(strcmp(UniqueDirectionLabels, oldLabels{i}));
            UniqueDirectionLabels{ix} = newLabels{i};
        end
    end
    % Create a vector that contains the Gaussian stimulation
    
    FWHM=3.8; % width of the stimulation Gaussian in seconds
    Sigma=FWHM/2.355; % Sigma of the Gaussian in seconds
    SigmaSamp=Sigma*params.sampling_frequency; % Sigma now in units of data samples
    Alpha=(params.full_trial_length*params.sampling_frequency)/(SigmaSamp*2);
    GaussModel=gausswin(params.full_trial_length*params.sampling_frequency,Alpha);
    
    GaussModel=fshift(GaussModel,round(params.full_trial_length*params.sampling_frequency/4));
    mkdir(fullfile(resultsPath, char(Subjects(SubjectID))));
    cd(fullfile(resultsPath, char(Subjects(SubjectID))));
    fprintf(['\n\t*** Averaging ' num2str(length(UniqueFreqs)*length(UniqueDirections)) ' crossings of frequency and direction: ']); % Update user
    
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
                theMean = mean(TimeSeriesMatrix(1:100, :));
                theMeanMatrix(Indices) = theMean;
                TimeSeriesMatrix = ((TimeSeriesMatrix-repmat(theMean, size(TimeSeriesMatrix, 1), 1))./repmat(theMean, size(TimeSeriesMatrix, 1), 1));
                theRemoveIdx = find(sum(abs(TimeSeriesMatrix) > 0.5));
                theRemoveIdxReverse = find(~(sum(abs(TimeSeriesMatrix) > 0.5)));
                TimeSeriesMatrix(:, theRemoveIdx) = [];
                theMean(:, theRemoveIdx) = [];
                theValid(Indices(theRemoveIdxReverse)) = 1;
                for m = 1:size(TimeSeriesMatrix, 2)
                    iy = TimeSeriesMatrix(:, m);
                    [iy, indx]=SpikeRemover(iy,params.spike_remover_params);
                    
                    removePoints = find(abs(iy)>params.BadPercentChangeThreshold);
                    iy(removePoints)=NaN;
                    TimeSeriesMatrix(:, m) = iy;
                end
            end
            MeanMatrix{f,d,:} = theMean;
            AvgTimeSeries(f,d,:)=nanmean(TimeSeriesMatrix,2);
            SEMTimeSeries(f,d,:)=nanstd(TimeSeriesMatrix, [], 2)/sqrt(size(TimeSeriesMatrix,2 ));
            TimeSeriesMatrixStore{f, d} = TimeSeriesMatrix;
            if isempty(strfind(UniqueDirectionLabels{d}, 'Background'));
                if params.SaveDataFlag
                    csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '_TimeSeries.csv'], TimeSeriesMatrix);
                end
            end
        end % if indices is not length zero
    end % for number of unique directions
    for j=0:log10(IterationCount-1)
        fprintf('\b'); % delete previous counter display
    end
    fprintf('done.'); % notify user we are done the loop
    
    timeVector = 0:1/params.sampling_frequency:params.final_trial_length-1/params.sampling_frequency;
    
    % Plot the time series
    for d=1:length(UniqueDirections)
        figure;
        hold on;
        %plot(timeVector(1:400), squeeze(AvgTimeSeries(:, d, 1:400)), '-k');
        plot([5+params.StepDurSecs 5+params.StepDurSecs], [-0.4 0.4], 'r'); hold on;
        shadedErrorBar(timeVector(1:600), squeeze(AvgTimeSeries(:, d, 1:600)), squeeze(SEMTimeSeries(:, d, 1:600)));
        plot([timeVector(1) timeVector(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
        pbaspect([1 1 1]);
        title({strrep(UniqueDirectionLabels{d}, '_', ' ') char(Subjects(SubjectID))});
        plot([5 5], [-0.4 0.4], 'r');
        ylim([-0.4 0.4]);
        M = [timeVector(1:600)' squeeze(AvgTimeSeries(:, d, 1:600))];
        if params.SaveDataFlag
            % Save out mean pupil size
            csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '.csv'], M);
            csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '_Mean.csv'], MeanMatrix{1, d,:});
        end
        if params.SavePlotFlag
            % Save out the plot
            set(gca, 'TickDir', 'out');
            set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
            set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
            saveas(gcf, [char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '.pdf'], 'pdf');
            close(gcf);
        end
    end
    fprintf('\n');
end % main
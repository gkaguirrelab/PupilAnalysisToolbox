function ReturnData = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath)
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
    clear invalid;
    
    % Prepare plot figures for this subject
    if (params.TrialInspectorFlag==1)
        figTrialInspector=figure;
        pause on;
    end
    
    % Create the results path
    fullResultsPath = fullfile(resultsPath, char(Subjects(SubjectID)));
    if ~exist(fullResultsPath, 'dir');
        mkdir(fullResultsPath);
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
            [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempTrialContrasts,TempAlphaSpacing,TempDateTime]=...
                PupilAnalysisToolbox_LoadData(basePath, char(Protocols(p)),char(Subjects(SubjectID)),sess);
            
            if (not(isempty(TempData)))
                fprintf('.'); % Notify user
                
                % Now remove those trials that have insufficient time points
                NumTrials = length(TempData);
                
                % Now concatenate TempData to the full Data array
                if (FirstGoodSessionFlag==1)
                    Data = TempData;
                    for trial = 1:NumTrials
                        Data(trial).frequencies = TempTrialFrequencies;
                        Data(trial).directions = TempTrialDirections;
                        Data(trial).phases = TempTrialPhases;
                        Data(trial).contrasts = TempTrialContrasts;
                    end
                    FirstGoodSessionFlag = 0;
                else
                    Data(end+1:end+NumTrials) = TempData(1:end);
                    Data(end+1:end+NumTrials).frequencies = TempTrialFrequencies(1:end);
                    Data(end+1:end+NumTrials).directions = TempTrialDirections(1:end);
                    Data(end+1:end+NumTrials).phases = TempTrialPhases(1:end);
                    Data(end+1:end+NumTrials).contrasts = TempTrialContrasts(1:end);
                end
            end % if file exists
        end % loop over directions
    end % loop over sessions
    
    fprintf('done.'); % Notify user
    
    % Smooth, interpolate, and de-spike the data
    fprintf(['\n\t*** Smoothing and interpolating ' num2str(length(Data)) ' trials...  ']); % Update user
    
    % Iterate over all the data
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
        if ~isempty(Data(trial).time)
            iy=SGolaySmooth(Data(trial).time, Data(trial).diameter,...
                params.sgolay_span, params.sgolay_polynomial, params.sampling_frequency,...
                params.full_trial_length);
            
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
            
            if (params.TrialInspectorFlag==1)
                figure(1);
                plot(Data(trial).time,Data(trial).diameter, '-r');
                figure(2);
                plot(iy, '-k');
                pbaspect([1 1 1]);
                pause;
                
            end
            
            % Align to zero stimulus phase.
            AmountToShift= TrialPhases(trial)*...  % Length of stimulus cycle in seconds
                params.sampling_frequency*...;  % number of data samples per second
                (1);    % we need to phase advance
            orig_iy = iy;
            tmp = NaN*ones(size(orig_iy));
            tmp(1:length((AmountToShift+1):length(iy))) = orig_iy((AmountToShift+1):length(iy));
        else
            tmp = NaN;
        end
        
        TimeSeries(trial, :) = tmp;
    end
    
    trialCounter = trial+1;
    
    % Trim the data variables to account for any skipped trials due to
    % exceeding the allowable proportion of NaNs
    TimeSeries = TimeSeries(1:trialCounter-1,:);
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
    % Clean up the direction labels
    UniqueDirections=unique(TrialDirections);
    for d = 1:length(UniqueDirections)
        tmp0 = allwords(UniqueDirections{d}, '-');
        UniqueDirectionLabels{d} = [tmp0{2} '-' tmp0{3}];
    end
    if (params.RelabelDirections==1)
        for i = 1:length(oldLabels)
            ix = strcmp(UniqueDirectionLabels, oldLabels{i});
            UniqueDirectionLabels{ix} = newLabels{i};
        end
    end
    
    fprintf(['\n\t*** Averaging ' num2str(length(UniqueFreqs)*length(UniqueDirections)) ' crossings of frequency and direction.']); % Update user
    dd = 0;
    
    timeVector = 0:1/params.sampling_frequency:params.final_trial_length-1/params.sampling_frequency;
    
    % Set up saving of two data quality files. One gives the summary
    % information regarding the number of rejected trials in each category
    % of trial type. The second provides a vector of which trials in the
    % order that they were acquired were marked as bad.
    dataQualityOutFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_DataQuality.csv']);
    fid_quality = fopen(dataQualityOutFile, 'w');
    fprintf(fid_quality, 'Direction,Rejected trials,Total trials,Percentage rejected\n');
    
    badTrialVectorOutFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_BadTrialVector.csv']);
    fid_badTrialVector = fopen(badTrialVectorOutFile, 'w');
    
    rejectNCumulative = 0;
    totalNCumulative = 0;
    for d=1:length(UniqueDirections)
        % Define the outfile
        timeSeriesOutFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '_TimeSeries.csv']);
        csvwrite(timeSeriesOutFile, TimeSeriesMatrix);
        
        IterationCount=(f-1)*length(UniqueDirections) + d;
        Indices = find(strcmp(TrialDirections,UniqueDirections(d)));
        
        dd = dd+1;
        if (not(isempty(Indices)))
            % Create the average time series. This is done by first
            %  assembling a matrix across time-series, and then using
            %  the nanmean to handle NaN points
            for i=1:length(Indices)
                if (i==1)
                    TimeSeriesMatrixOrig=[TimeSeries(Indices(i),:)'];
                else
                    TimeSeriesMatrixOrig=[TimeSeriesMatrixOrig TimeSeries(Indices(i),:)'];
                end % if the first timeseries
            end % for number of indicies
            theMean = mean(TimeSeriesMatrixOrig(1:100, :));
            theMeanMatrix(Indices) = theMean;
            TimeSeriesMatrix = ((TimeSeriesMatrixOrig-repmat(theMean, size(TimeSeriesMatrixOrig, 1), 1))./repmat(theMean, size(TimeSeriesMatrixOrig, 1), 1));
            TimeSeriesMatrixOrigNonClean = TimeSeriesMatrixOrig;
            for m = 1:size(TimeSeriesMatrix, 2)
                iy = TimeSeriesMatrix(:, m);
                %[iy, indx] = SpikeRemover(iy,params.spike_remover_params);
                
                removePoints = find(abs(iy) > params.BadPercentChangeThreshold);
                iy(removePoints) = NaN;
                velocityTrace = diff(smooth(iy, params.VelocitySmoothingParam));
                % Go through the velocity signal and try to find the
                % onset of a blink event
                foundPeaks = [];
                
                v = 0;
                while v < length(velocityTrace)
                    v = v+1;
                    if velocityTrace(v) < params.VelocityOnsetThreshold
                        onsetCandidatePeak = v;
                        % Search ahead in the window to find if we
                        % exceed the offset threshold also
                        if v+params.VelocitySearchWindowSize > length(velocityTrace)
                            searchIndices = v:length(velocityTrace);
                        else
                            searchIndices = v:v+params.VelocitySearchWindowSize;
                        end
                        for vv = searchIndices
                            if velocityTrace(vv) > params.VelocityOffsetThreshold
                                offsetCandidatePeak = vv;
                                v = offsetCandidatePeak+params.VelocityMarginWindowSize;
                                
                                if (onsetCandidatePeak-params.VelocityMarginWindowSize) < 1
                                    startIdx = 1;
                                else
                                    startIdx = (onsetCandidatePeak-params.VelocityMarginWindowSize);
                                end
                                if (offsetCandidatePeak+params.VelocityMarginWindowSize) > length(velocityTrace)
                                    endIdx = length(velocityTrace);
                                else
                                    endIdx = (offsetCandidatePeak+params.VelocityMarginWindowSize);
                                end
                                foundPeaks = [foundPeaks startIdx:endIdx];
                                break;
                            end
                        end
                        
                    end
                    iy(foundPeaks) = NaN;
                    TimeSeriesMatrix(:, m) = iy;
                    TimeSeriesMatrixOrig(foundPeaks, m) = NaN;
                    TimeSeriesMatrixOrig(removePoints, m) = NaN;
                end
                
                v = 0;
                while v < length(velocityTrace)
                    v = v+1;
                    if velocityTrace(v) > params.VelocityOffsetThreshold
                        onsetCandidatePeak = v;
                        % Search ahead in the window to find if we
                        % exceed the offset threshold also
                        if v+params.VelocitySearchWindowSize > length(velocityTrace)
                            searchIndices = v:length(velocityTrace);
                        else
                            searchIndices = v:v+params.VelocitySearchWindowSize;
                        end
                        for vv = searchIndices
                            if velocityTrace(vv) < params.VelocityOnsetThreshold
                                offsetCandidatePeak = vv;
                                v = offsetCandidatePeak+params.VelocityMarginWindowSize;
                                if (onsetCandidatePeak-params.VelocityMarginWindowSize) < 1
                                    startIdx = 1;
                                else
                                    startIdx = (onsetCandidatePeak-params.VelocityMarginWindowSize);
                                end
                                if (offsetCandidatePeak+params.VelocityMarginWindowSize) > length(velocityTrace)
                                    endIdx = length(velocityTrace);
                                else
                                    endIdx = (offsetCandidatePeak+params.VelocityMarginWindowSize);
                                end
                                foundPeaks = [foundPeaks startIdx:endIdx];
                                break;
                            end
                        end
                    end
                    iy(foundPeaks) = NaN;
                    TimeSeriesMatrix(:, m) = iy;
                    TimeSeriesMatrixOrig(foundPeaks, m) = NaN;
                    TimeSeriesMatrixOrig(removePoints, m) = NaN;
                end
                
                if params.TrialInspectorFlag
                    figure(figTrialInspector);
                    hold off;
                    plot(TimeSeriesMatrixOrigNonClean(:, m), '-r');  hold on;
                    plot(TimeSeriesMatrixOrig(:, m), '-k');
                    pbaspect([1 1 1]);
                    pause;
                end
                
            end
            
            % Remove points beyond what is in minimum_length_trial
            TimeSeriesMatrixOrig = TimeSeriesMatrixOrig(1:params.sampling_frequency*params.minimum_length_trial, :);
            
            % Make sure we're within the bounds for rejecting an entire
            % trial
            totalN(d) = size(TimeSeriesMatrixOrig, 2);
            rejectN(d) = 0;
            for m = 1:size(TimeSeriesMatrixOrig, 2)
                if sum(isnan(TimeSeriesMatrixOrig(:, m)))/length(TimeSeriesMatrixOrig(:, m)) > params.BadNaNThreshold;
                    invalid(m) = 1;
                    rejectN(d) = rejectN(d)+1;
                    fprintf(fid_badTrialVector,'1');
                else
                    invalid(m) = 0;
                    fprintf(fid_badTrialVector,'0');
                end
                
                if m ~= size(TimeSeriesMatrixOrig, 2)
                    fprintf(fid_badTrialVector,',');
                end
            end
            rejectNCumulative = rejectNCumulative+rejectN(d);
            totalNCumulative = totalNCumulative+totalN(d);
            
            
            invalid = logical(invalid);
            TimeSeriesMatrixOrig(:, invalid) = [];
            
            % Save and report the data quality information
            fprintf(fid_quality, '%s,%g,%g,%.2f\n', UniqueDirectionLabels{d}, rejectN(d), totalN(d), (rejectN(d)/totalN(d)));
            fprintf('\n\t-> %s: rejecting %g of %g (%.2f)', UniqueDirectionLabels{d}, rejectN(d), totalN(d), (rejectN(d)/totalN(d)));
            
            
            % Now, with the blinks etc. removed, do normalize by mean
            % again
            theMean = mean(TimeSeriesMatrixOrig(1:100, :));
            TimeSeriesMatrix = ((TimeSeriesMatrixOrig-repmat(theMean, size(TimeSeriesMatrixOrig, 1), 1))./repmat(theMean, size(TimeSeriesMatrixOrig, 1), 1));
        end
        
        MeanMatrix{f,d,:} = theMean;
        AvgTimeSeries(f,d,:) = nanmean(TimeSeriesMatrix,2);
        SEMTimeSeries(f,d,:) = nanstd(TimeSeriesMatrix, [], 2)/sqrt(size(TimeSeriesMatrix,2 ));
        TimeSeriesMatrixStore{f, d} = TimeSeriesMatrix;
        
        ReturnData(SubjectID, dd).label = UniqueDirectionLabels{d};
        ReturnData(SubjectID, dd).t = timeVector(1:length(AvgTimeSeries));
        ReturnData(SubjectID, dd).Mean = theMean;
        ReturnData(SubjectID, dd).TimeSeries = TimeSeriesMatrix;
        ReturnData(SubjectID, dd).AvgTimeSeries = nanmean(TimeSeriesMatrix,2);
        ReturnData(SubjectID, dd).SEMTimeSeries = nanstd(TimeSeriesMatrix, [], 2)/sqrt(size(TimeSeriesMatrix,2 ));
        
    end % if indices is not length zero
    fprintf(fid_quality, 'Total,%g,%g,%.2f\n', rejectNCumulative, totalNCumulative, (rejectNCumulative/totalNCumulative));
    
end % for number of unique directions
fclose(fid_quality);
fclose(fid_badTrialVector);

xLim = params.xLim;
yLim = params.yLim;
% Plot the time series
for d=1:length(UniqueDirections)
    if ~strcmp(UniqueDirections{d}, 'Background') % Only do this for non-background data
        figure;
        hold on;
        %plot(timeVector(1:400), squeeze(AvgTimeSeries(:, d, 1:400)), '-k');
        hold on;
        shadedErrorBar(timeVector(1:length(AvgTimeSeries)), 100*squeeze(AvgTimeSeries(:, d, :)), 100*squeeze(SEMTimeSeries(:, d, :)));
        plot([timeVector(1) timeVector(length(AvgTimeSeries))], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
        pbaspect([1 1 1]);
        title({strrep(UniqueDirectionLabels{d}, '_', ' ') strrep(char(Subjects(SubjectID)), '_', ' ')});
        plot([5 5], 100*[-yLim yLim], 'r');
        ylim(100*[-yLim yLim]);
        xlim([0 xLim]);
        M = [timeVector(1:500)' squeeze(AvgTimeSeries(:, d, 1:500))];
        xlabel('Time [s]');
        ylabel('\Delta%');
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
            saveas(gcf, [char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '.png'], 'png');
            close(gcf);
        end
    end
    fprintf('\n');
end % main

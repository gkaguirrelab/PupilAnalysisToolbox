function ReturnData = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocols, newLabels, oldLabels, basePath, resultsPath)
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
    clear rData;
    clear pData;
    
    % Prepare plot figures for this subject
    if (params.TrialInspectorFlag == 1)
        figTrialInspector = figure;
    end
    
    % Create the results path
    fullResultsPath = fullfile(resultsPath, char(Subjects(SubjectID)), char(Dates(SubjectID)));
    if ~exist(fullResultsPath, 'dir');
        mkdir(fullResultsPath);
    end
    
    % Display some information
    fprintf('\n>>> Processing subject <strong>%s</strong>', char(Subjects(SubjectID)));
    fprintf('\n\t*** Loading and concatenating sessions for subject <strong>%s</strong>', char(Subjects(SubjectID)));
    
    % Loop over sessions and directions for a subject and concatenate the
    % data. The number of sessions is currently hardcoded a 3, as we have
    % not implemented a mechanism for detecting how many files across
    % sessions are present.
    
    FirstGoodSessionFlag = 1;
    for sess = 1:5
        for p = 1:length(Protocols)
            % Load data set
            [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempTrialContrasts,TempAlphaSpacing,TempDateTime]=...
                PupilAnalysisToolbox_LoadData(basePath, char(Protocols(p)),char(Subjects(SubjectID)),char(Dates(SubjectID)),sess);
            
            if (not(isempty(TempData)))
                fprintf('.'); % Notify user
                
                % In prior analyses, this is the point at which trials
                % would be discarded for not having a sufficient number of
                % data points. This step is no longer performed here.
                % Instead, we will simply set NumTrials to be equal to the
                % total number of trials stored in TempData
                
                NumTrials = length(TempData);
                
                % Now concatenate TempData to the full Data array
                if (FirstGoodSessionFlag==1)
                    rData = TempData;
                    FirstGoodSessionFlag = 0;
                    for trial = 1:NumTrials
                        % Assign the direction label
                        tmp0 = allwords(TempTrialDirections{trial}, '-');
                        rData(trial).direction = [tmp0{2} '-' tmp0{3}];
                    end
                else
                    
                    % On August 9, 2016, gka switched the row-column order
                    % of the size function for the check on rData. At some
                    % point prior to this date, the row-column order of the
                    % rData variable was switched, causing the aggregation
                    % of data from across sessions to not work. The change
                    % of the code from "size(rData,1)" to "size(rData,2)"
                    % fixed this problem. Future Geoff, you may reach this
                    % comment when you are trying to analyze old pupil data
                    % and discover that the aggregation of data from across
                    % sessions is not working.
                    
                    endIdx = size(rData, 2);
                    rData(endIdx+1:endIdx+NumTrials) = TempData(1:end);
                    for trial = 1:NumTrials
                        % Assign the direction label
                        tmp0 = allwords(TempTrialDirections{trial}, '-');
                        rData(endIdx+trial).direction = [tmp0{2} '-' tmp0{3}];
                    end
                end
            end % if file exists
        end % loop over directions
    end % loop over sessions
    
    % Clean up the data struct
    pData = rData;
    pData = rmfield(pData, 'time_inter');
    pData = rmfield(pData, 'modulationMode');
    pData = rmfield(pData, 'frequencyCarrier');
    pData = rmfield(pData, 'frequencyEnvelope');
    pData = rmfield(pData, 'phaseCarrier');
    pData = rmfield(pData, 'phaseEnvelope');
    pData = rmfield(pData, 'rawMmPositions');
    try
        pData = rmfield(pData, 'rawFickPositions');
    end
    try
        pData = rmfield(pData, 'rawFickPsotions');
    end
    pData = rmfield(pData, 'rawHelmholtzPositions');
    [pData.timeMSecsRaw] = pData.rawTimeStamps;
    pData = rmfield(pData, 'time');
    pData = rmfield(pData, 'rawTimeStamps');
    [pData.pupilDiameterMmRaw] = pData.rawPupilDiameter;
    [pData.pupilDiameterMm] = pData.rawPupilDiameter;
    pData = rmfield(pData, 'diameter');
    pData = rmfield(pData, 'rawPupilDiameter');
        
    %% Remove the first 5 trials as they are background adaptation
    
    % This needs a parameter setting to adjust depending on the
    % implementation.
    
%    pData(1:5) = [];
%    rData(1:5) = [];
    NumTrials = length(pData);
    
    fprintf('- Done.\n'); % Notify user

    % Interpolate, and de-spike the data
    fprintf('\t*** Smoothing and interpolating <strong>%g</strong> trials...  ', NumTrials); % Update user
    
    % Iterate over all the data
    for trial = 1:NumTrials
        
        % Trial counter feedback to user
        for j = 0:log10(trial-1)
            fprintf('\b');          % delete previous counter display
        end
        fprintf('%d', trial);       % update progress trial number
        
        % Shift the time stamps in the data array back by the
        % params.StimOnsetDelay value. This accounts for the delay in the onset of
        % the stimulus modulation relative to the time stamps on the pupil
        % size measurements. This will result in some time stamps having a
        % negative value, but the subsequent SGolaySmooth step will handle
        % this and return a vector that is just from time point zero
        % onward.
        % Replace all 0 with NaN
        
        pData(trial).pupilDiameterMm(pData(trial).pupilDiameterMm == 0) = NaN;
        pData(trial).pupilDiameterMm = pData(trial).pupilDiameterMm';
        pData(trial).timeMSecs = pData(trial).timeMSecsRaw'-pData(trial).timeMSecsRaw(1);
        pData(trial).timeMSecs = pData(trial).timeMSecs-params.StimOnsetDelay;
        
        % Clip and phase-correct the data.
        % Check  that the trial is not empty, or contains a single
        % value. If there is trial data present, proceed.
        if ~isempty(pData(trial).timeMSecs) && (~(numel(pData(trial).timeMSecs) == 1))
            
            % the variable iy holds the pupil diameter data that will then
            % be manipulated.
            
            iy = pData(trial).pupilDiameterMm;
            
            % Clip the vector to correspond to the length of the full
            % data for the trial: params.full_trial_length
            iy = iy(1:params.full_trial_length*params.sampling_frequency);
            
            % Temporarily crop to params.final_trial_length to allow corrections of
            % phase offsets in the stimuli.
            % Check to make sure that we are not asking for a final trial
            % length longer than the full trial length
            if (params.final_trial_length > params.full_trial_length)
                error('params.final_trial_length must be less than or equal to params.full_trial_length')
            end
            
            % Align to zero stimulus phase.
            AmountToShift = pData(trial).phaseRandSec*...  % Length of stimulus cycle in seconds
                params.sampling_frequency*...;  % number of data samples per second
                (1);    % we need to phase advance
            orig_iy = iy;
            tmp = NaN*ones(size(orig_iy));
            tmp(1:length((AmountToShift+1):length(iy))) = orig_iy((AmountToShift+1):length(iy));
            
            pData(trial).timeSecs = linspace(0, params.minimum_trial_length-(1/params.sampling_frequency), params.minimum_trial_length*params.sampling_frequency);
            pData(trial).timeMSecs = 1000*pData(trial).timeSecs;
            pData(trial).pupilDiameterMm = tmp(1:length(pData(trial).timeSecs));
        else % The trial is full of NaNs or has a single value. Manufacture
            % a place-holder trial of the appropriate length that is all NaNs
            pData(trial).timeSecs = linspace(0, params.minimum_trial_length-(1/params.sampling_frequency), params.minimum_trial_length*params.sampling_frequency);
            pData(trial).timeMSecs = 1000*pData(trial).timeSecs;
            pData(trial).pupilDiameterMm = NaN*ones(size(pData(trial).timeSecs));
        end
    end
    fprintf('. - Done.\n'); % notify user we are done the loop
    
    %% Data quality: Spike removing, making sure we have enough data points, ...
    fprintf('\t*** De-spiking and mean-centering <strong>%g</strong> trials... ', NumTrials); % Update user
    
    % Window to be used to calculate where to set zero on the y-axis
    meanCenterWindowIdx = find((pData(trial).timeSecs >= params.meanCenterWindow(1)) & (pData(trial).timeSecs <= params.meanCenterWindow(2)));
    
    for trial = 1:NumTrials
        % Trial counter information for the user
        for j=0:log10(trial-1)
            fprintf('\b');          % delete previous counter display
        end
        fprintf('%d', trial);       % update progress trial number
        
        % 1. Mean-center the data
        pData(trial).meanBaseline = nanmean(pData(trial).pupilDiameterMm(meanCenterWindowIdx));
        pData(trial).pupilDiameterMmMeanCentered = ((pData(trial).pupilDiameterMm-repmat(pData(trial).meanBaseline, ...
            size(pData(trial).pupilDiameterMm, 1), 1))./repmat(pData(trial).meanBaseline, ...
            size(pData(trial).pupilDiameterMm, 1), 1));
        
        % 2. Get the indices from the spikes
        [~, pData(trial).removePoints] = PupilAnalysisToolbox_SpikeRemover(pData(trial).pupilDiameterMmMeanCentered, params);
        
        % 3. Remove the spikes from the original data (set them to NaN)
        pData(trial).pupilDiameterMmDespiked = pData(trial).pupilDiameterMm;
        pData(trial).pupilDiameterMmDespiked(pData(trial).removePoints) = NaN;
        
        % 4. Interpolate
        theNans = isnan(pData(trial).pupilDiameterMmDespiked);
        if sum(~theNans) > 1
            x = pData(trial).pupilDiameterMmDespiked;
            x(theNans) = interp1(pData(trial).timeSecs(~theNans), pData(trial).pupilDiameterMmDespiked(~theNans), ...
                pData(trial).timeSecs(theNans), 'linear');
            pData(trial).pupilDiameterMmDespiked = x;
        end
        
        % 5. Mean-center again
        pData(trial).meanBaseline = nanmean(pData(trial).pupilDiameterMmDespiked(meanCenterWindowIdx));
        pData(trial).pupilDiameterMmMeanCentered = ((pData(trial).pupilDiameterMmDespiked-repmat(pData(trial).meanBaseline, ...
            size(pData(trial).pupilDiameterMm, 1), 1)) ./ repmat(pData(trial).meanBaseline, ...
            size(pData(trial).pupilDiameterMm, 1), 1));
        
        % 6. Relabel the direction name
        labelIdx = find(ismember(lower(oldLabels), lower(pData(trial).direction)));
        pData(trial).direction = newLabels{labelIdx};
        
        % 7. Determine the proportion of missing data points
        pData(trial).propMissingData = sum(theNans)/length(pData(trial).pupilDiameterMmDespiked);
        if pData(trial).propMissingData > params.BadNaNThreshold
            pData(trial).dataQualityPass = 0;
        else
            pData(trial).dataQualityPass = 1;
        end
        
        % 8. (optional) Display the trial and the removed points
        if (params.TrialInspectorFlag == 1)
            figure(figTrialInspector); hold off;
            plot(pData(trial).timeSecs, pData(trial).pupilDiameterMm, '-r'); hold on;
            plot(pData(trial).timeSecs, pData(trial).pupilDiameterMmDespiked, '-k');
            xlabel('Time [Secs]'); ylabel('Pupil diameter [mm]');
            xlim([0 params.final_trial_length]); ylim([2 9]);
            pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
            if pData(trial).dataQualityPass == 1
                title({['Trial ' num2str(trial)] 'PASS' num2str(100*pData(trial).propMissingData)});
            else
                title({['Trial ' num2str(trial)] 'REJECT' num2str(100*pData(trial).propMissingData)});
            end
            pause;
        end
    end
    fprintf('. - Done.\n'); % notify user we are done the loop
    
    %% Data quality
    % 1. Trial-by-trial information
    badTrialVectorOutFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_BadTrialVector.csv']);
    fid = fopen(badTrialVectorOutFile, 'w');
    fprintf(fid, 'Trial #,Direction,Prop. missing data,Pass?\n');
    for trial = 1:NumTrials
        fprintf(fid, '%g,%s,%.3f,%g\n', trial, pData(trial).direction, pData(trial).propMissingData, pData(trial).dataQualityPass);
    end
    fclose(fid);
    
    % 2. Information by direction
    dataQualityOutFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_DataQuality.csv']);
    fid = fopen(dataQualityOutFile, 'w');
    fprintf(fid, 'Direction,Rejected trials,Total trials,Percentage rejected\n');
    for dd = 1:length(newLabels)
        ind = find(ismember({pData.direction}, newLabels{dd}));
        rejectedTrials = sum(~[pData(ind).dataQualityPass]);
        totalTrials = length(ind);
        percentageRejected = rejectedTrials / totalTrials;
        fprintf(fid, '%s,%g,%g,%.3f\n', newLabels{dd}, rejectedTrials, totalTrials, percentageRejected);
    end
    fclose(fid);
    
    %% Calculate averages
    for dd = 1:length(newLabels)
        
        % Identify good trials for this modulation, check if there are any
        ind = find(ismember({pData.direction}, newLabels{dd}) & [pData.dataQualityPass]);
        if ~(isempty(ind))
            
            % 1. Assemble the data
            ReturnData(SubjectID, dd).label = newLabels{dd};
            ReturnData(SubjectID, dd).timeSecs = pData(ind(1)).timeSecs-params.PulseOnsetSecs';
            ReturnData(SubjectID, dd).TimeSeries = cell2mat({pData(ind).pupilDiameterMmMeanCentered}')';
            ReturnData(SubjectID, dd).Mean = [pData(ind).meanBaseline];
            ReturnData(SubjectID, dd).AvgTimeSeries = nanmean(ReturnData(SubjectID, dd).TimeSeries, 2);
            ReturnData(SubjectID, dd).SEMTimeSeries = nanstd(ReturnData(SubjectID, dd).TimeSeries, [], 2) / sqrt(size(ReturnData(SubjectID, dd).TimeSeries, 2));
            
            % 2. Plot the data
            plot([-params.PulseOnsetSecs params.xLim], [0 0], '-', 'Color', [0.3 0.3 0.3]); hold on;
            shadedErrorBar(ReturnData(SubjectID, dd).timeSecs, ReturnData(SubjectID, dd).AvgTimeSeries, ...
                ReturnData(SubjectID, dd).SEMTimeSeries);
            plot([0 params.PulseDurationSecs], [0.1 0.1], '-r', 'LineWidth', 2);
            pbaspect([1 1 1]);
            xlim([-params.PulseOnsetSecs params.xLim]); ylim([-params.yLim params.yLim]);
            xlabel('Time [Secs]'); ylabel('Pupil diameter [change]');
            pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
            title({newLabels{dd} ['Mean\pm1SEM (n = ' num2str(length(ind)) ' trials)']});
            
            % 3. Save the plot
            if params.SavePlotFlag
                % Save out the plot
                set(gca, 'TickDir', 'out');
                set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 4 and height 4.
                set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 4 and height 4.
                outFile1 = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_' newLabels{dd} '.pdf']);
                saveas(gcf, outFile1, 'pdf');
                outFile2 = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_' newLabels{dd} '.png']);
                saveas(gcf, outFile2, 'png');
                close(gcf);
            end
            
            % 4. (optional) Save the data
            if params.SaveDataFlag
                % Save out mean pupil size
                outFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_' newLabels{dd} '_TimeSeries.csv']);
                csvwrite(outFile, ReturnData(SubjectID, dd).TimeSeries);
                outFile = fullfile(fullResultsPath, [char(Subjects(SubjectID)) '_PupilPulseData_' newLabels{dd} '_Mean.csv']);
                csvwrite(outFile, ReturnData(SubjectID, dd).Mean);
            end
            
        else % There were no good trials, return NaNs
            ReturnData(SubjectID, dd).label = newLabels{dd};
            ReturnData(SubjectID, dd).timeSecs = linspace(0, params.minimum_trial_length-(1/params.sampling_frequency), params.minimum_trial_length*params.sampling_frequency)';
            ReturnData(SubjectID, dd).TimeSeries = NaN*ones(size(pData(trial).timeSecs))';
            ReturnData(SubjectID, dd).Mean = NaN;
            ReturnData(SubjectID, dd).AvgTimeSeries = NaN*ones(size(pData(trial).timeSecs))';
            ReturnData(SubjectID, dd).SEMTimeSeries = NaN*ones(size(pData(trial).timeSecs))';
        end
    end
    
end
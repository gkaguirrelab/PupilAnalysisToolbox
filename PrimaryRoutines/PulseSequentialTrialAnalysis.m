function [TimeSeries]=PulseSequentialTrialAnalysis(...
    sampling_frequency,...
    full_trial_length,...
    final_trial_length,...
    minimum_length_trial,...
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
    StimOnsetDelay);



% Add the OLPupilDiameter folder to the path. Different conventions.
if isdir('/Users/Shared/Matlab/Experiments/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OLPupilDiameter';
elseif isdir('/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter';
end

ResultsFullPath = fullfile(basePath, 'analysis', 'results');

% Calculate the  length of our data vectors

datalen = full_trial_length*sampling_frequency;

% Set up some figure windows

maxSubplots = length(Subjects);
figMeanPupilByDirection = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RETRIEVE CONTRASTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to obtain the contrasts that each of the photoreceptors gets.
% This is done with a call to OLPDCalculateContrast(...).
configPath = fullfile(basePath, 'code', 'config');

% Load the config files. We hard code 0 in there. This is the index for
% which we want to calculate contrasts. That is 90 deg phase angle when we
% have a stimulus titration with 200 steps between 0 and 2*pi.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ROUTINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now loop across subjects and perform the analysis

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
    
    
    
    
    % Prepeare plot figures for this subject
    
    if (TrialInspectorFlag==1)
        figTrialInspector=figure;
        pause on;
    end
    
    figSparkLinePlots(SubjectID) = figure;
    
    SubjectAnon = ['S' num2str(SubjectID)];
    
    fprintf(['\nLoading and concatenating sessions for subject '...
        char(Subjects(SubjectID))]); % Notify user
    
    % Loop over sessions and directions for a subject and concatenate the data
    
    FirstGoodSessionFlag=1;
    
    for sess=1:3
        for p=1:length(Protocols)
            
            % Load data set
            [TempData,TempTrialFrequencies,TempTrialPhases,TempTrialDirections,TempAlphaSpacing,TempDateTime]=...
                DataLoaderFlickerSensitivity(char(Protocols(p)),char(Subjects(SubjectID)),sess);
            
            if (not(isempty(TempData)))
                
                fprintf('.'); % Notify user
                
                % Keep AlphaSpacing
                
                AlphaSpacing=TempAlphaSpacing;
                
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
                    if ((TempData(trial).time(end)/1000) < minimum_length_trial)
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
    
    fprintf('done\n'); % Notify user
    
    % Define some variables to hold the processed data
    
    %     TimeSeries=zeros(length(Data),datalen);
    %     MeanPupilByTrial=zeros(1,length(Data));
    %     AmplitudesByTrial=zeros(1,length(Data));
    %     PhasesByTrial=zeros(1,length(Data));
    
    % Smooth, interpolate, and de-spike the data
    %  Calculate and store mean pupil, amplitude, and phase
    
    fprintf(['\nSmoothing and interpolating ' num2str(length(Data)) ' trials:  ']); % Update user
    
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
        % phase value is the step on which the stimulus started. For the
        % Gaussian pulse stimuli, the phase shift is performed simply as a
        % shift in the specified mirror settings. The shift direction
        % is always positive, as the value in TrialPhases is either zero
        % (requiring no shift) or a positive integer.
        
        % One is subtracted from the TrialPhases value as a value of one
        % corresponds to the first mirror setting group, which is a phase
        % advancement of zero.
        
        AmountToShift= TrialPhases(trial)*...  % Length of stimulus cycle in seconds
            sampling_frequency*...;  % number of data samples per second
            (1);    % we need to phase advance
        
        % fshift implements a sinc shift, as we may have a non-integer shift
        orig_iy = iy;
        %iy=fshift(orig_iy,AmountToShift);
        tmp = NaN*ones(size(orig_iy));
        tmp(1:length((AmountToShift+1):length(iy))) = orig_iy((AmountToShift+1):length(iy));
        TimeSeries(trial, :) = tmp;
        % Store the mean pupil size.
        
        if (TrialInspectorFlag==1)
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
    
    if ~(TrialInspectorFlag==1)
        for j=0:log10(trial-1)
            fprintf('\b'); % delete previous counter display
        end
    end
    
    fprintf('done\n'); % notify user we are done the loop
    
    fprintf(['Found ' num2str(trialCounter-1) ' good trials.']);
    
    % Allocate variables to hold results
    
    % Create an average time-series for each unique frequency and direction
    % crossing. Calculate for these the best fit, amplitude and phase.
    % Also, calculate a noise level bootstrap for each crossing.
    
    UniqueFreqs=unique(TrialFrequencies);
    UniqueDirections=unique(TrialDirections);
    
    % Prepare some variables to hold the results
    
    
    % Create a vector that contains the Gaussian stimulation
    
    FWHM=3.8; % width of the stimulation Gaussian in seconds
    Sigma=FWHM/2.355; % Sigma of the Gaussian in seconds
    SigmaSamp=Sigma*sampling_frequency; % Sigma now in units of data samples
    Alpha=(full_trial_length*sampling_frequency)/(SigmaSamp*2);
    GaussModel=gausswin(full_trial_length*sampling_frequency,Alpha);
    
    % A hard-coded kludge. The Gaussian pulse was shifted back one-quarter
    % of the total duration of the trial
    UniqueDirectionLabels = UniqueDirections;
    if (RelabelDirections==1)
        for i = 1:length(oldLabels)
            ix = find(strcmp(UniqueDirectionLabels, oldLabels{i}));
            UniqueDirectionLabels{ix} = newLabels{i};
        end
    end
    GaussModel=fshift(GaussModel,round(full_trial_length*sampling_frequency/4));
    cd(ResultsFullPath)
    mkdir(ResultsDirName)
    cd(ResultsDirName);
    fprintf(['\nAveraging ' num2str(length(UniqueFreqs)*length(UniqueDirections)) ' crossings of frequency and direction: ']); % Update user
    
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
                    [iy, indx]=SpikeRemover(iy,spike_remover_params);
                    
                    removePoints = find(abs(iy)>BadPercentChangeThreshold);
                    iy(removePoints)=NaN;
                    TimeSeriesMatrix(:, m) = iy;
                end
                
                
            end
            
            MeanMatrix{f,d,:} = theMean;
            AvgTimeSeries(f,d,:)=nanmean(TimeSeriesMatrix,2);
            SEMTimeSeries(f,d,:)=nanstd(TimeSeriesMatrix, [], 2)/sqrt(size(TimeSeriesMatrix,2 ));
            TimeSeriesMatrixStore{f, d} = TimeSeriesMatrix;
            
            if isempty(strfind(UniqueDirectionLabels{d}, 'Background'));
                csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '_TimeSeries.csv'], TimeSeriesMatrix);
            end
        end % if indices is not length zero
    end % for number of unique directions
    
    for j=0:log10(IterationCount-1)
        fprintf('\b'); % delete previous counter display
    end
    fprintf('done\n'); % notify user we are done the loop
    
    
    
    timeVector = 0:1/sampling_frequency:final_trial_length-1/sampling_frequency;
    
    for d=2:length(UniqueDirections)
        subplot(2, ceil((length(UniqueDirections)-1)/2), d-1);
        hold on;
        %plot(timeVector(1:400), squeeze(AvgTimeSeries(:, d, 1:400)), '-k');
        shadedErrorBar(timeVector(1:600), squeeze(AvgTimeSeries(:, d, 1:600)), squeeze(SEMTimeSeries(:, d, 1:600)));
        plot([timeVector(1) timeVector(600)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
        pbaspect([1 1 1]);
        title(strrep(UniqueDirectionLabels{d}, '_', ' '));
        plot([5 5], [-0.4 0.4], 'r');
        plot([10 10], [-0.4 0.4], 'r'); hold on;
        ylim([-0.4 0.4]);
        
        M = [timeVector(1:600)' squeeze(AvgTimeSeries(:, d, 1:600))];
        
        csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '.csv'], M);
        
        
        % Save out mean pupil size
        csvwrite([char(Subjects(SubjectID)) '_PupilPulseData_' UniqueDirectionLabels{d} '_Mean.csv'], MeanMatrix{1, d,:});
    end
    
    set(gcf, 'PaperPosition', [0 0 12 6]); %Position plot at left hand corner with width 15 and height 6.
    set(gcf, 'PaperSize', [12 6]); %Set the paper to have width 15 and height 6.
    
    
    
    saveas(gcf, [char(Subjects(SubjectID)) '_PupilPulseData.pdf'], 'pdf');
    close(gcf);
    
    
    
    %
    %     shadedErrorBar(timeVector(1:400), squeeze(AvgTimeSeries(:, 3, 1:400)), squeeze(SEMTimeSeries(:, 3, 1:400)), {'Color', 'r'}); hold on;
    %     shadedErrorBar(timeVector(1:400), squeeze(AvgTimeSeries(:, 5, 1:400)), squeeze(SEMTimeSeries(:, 5, 1:400)), {'Color', 'b'});
    %         shadedErrorBar(timeVector(1:400), squeeze(AvgTimeSeries(:, d, 1:400)), squeeze(SEMTimeSeries(:, d, 1:400)));
    %     plot([timeVector(1) timeVector(400)], [0 0]', '-', 'Color', [0.2 0.2 0.2]);
    %                     plot([5 5], [-0.1 0.15], 'k');
    %             plot([10 10], [-0.15 0.15], 'k'); hold on;
    %                         pbaspect([1 1 1]);
    %                         ylim([-0.13 0.13]);
    %                 set(gcf, 'PaperPosition', [0 0 4 4]); %Position plot at left hand corner with width 15 and height 6.
    %     set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 15 and height 6.
    %
    %     saveas(gcf, 'PupilPulseDataLMSMel.pdf', 'pdf');
    %     close(gcf);
    %     fprintf('\nPlot results.\n'); % Notify user
    %
    %
    %     % Make a different label vector for plotting purposes
    %     UniqueDirectionLabels = UniqueDirections;
    %     if (RelabelDirections==1)
    %         for i = 1:length(oldLabels)
    %             ix = find(strcmp(UniqueDirectionLabels, oldLabels{i}));
    %             UniqueDirectionLabels{ix} = newLabels{i};
    %         end
    %     end
    %
    %     %% Create the sparkline plots
    %
    %     ErrorMultiplier=2; % Plots ±2SEM on time series sparklines
    %
    %     % first need to determine a constant Y-range across all plot cells
    %     for f=1:length(UniqueFreqs)
    %         for d=1:length(UniqueDirections)
    %             MaxYRange(f,d) = nanmax([ nanmax(squeeze(AvgTimeSeries(f,d,:))+ErrorMultiplier*squeeze(SEMTimeSeries(f,d,:))) nanmax(squeeze(AvgModelFits(f,d,:)))])  ;
    %             MinYRange(f,d) = nanmin([ nanmin(squeeze(AvgTimeSeries(f,d,:))-ErrorMultiplier*squeeze(SEMTimeSeries(f,d,:))) nanmin(squeeze(AvgModelFits(f,d,:)))])  ;
    %         end
    %     end
    %     yRange=[min(min(MinYRange)) max(max(MaxYRange))];
    %     yRange(1)=floor(yRange(1)*100)/100;
    %     yRange(2)=ceil(yRange(2)*100)/100;
    %
    %     % Now make the plots. One column for each direction
    %
    %     figure(figSparkLinePlots(SubjectID));
    %     SuperPlotTitle=['Subject ' char(SubjectAnon) '/Responses and fits [' num2str(yRange(1)) ',' num2str(yRange(2)) ']'];
    %     annotation('textbox', [0 0.9 1 0.1], ...
    %         'String', SuperPlotTitle, ...
    %         'EdgeColor', 'none', ...
    %         'HorizontalAlignment', 'center', ...
    %         'FontSize',12)
    %     if (and(length(UniqueDirections)==1,length(UniqueFreqs)==1))
    %         PlotTitle=char(UniqueDirectionLabels(1));
    %         RowLabels=num2str(UniqueFreqs);
    %         PlotOLPupilSparklines(AvgTimeSeries(:,d,:),AvgModelFits(:,d,:),SEMTimeSeries(:,d,:), ErrorMultiplier, yRange, figSparkLinePlots(SubjectID), RowLabels, PlotTitle, 1, 1);
    %     else % if we only have one crossing of freq and direction, else...
    %
    %         for d=1:length(UniqueDirections)
    %             PlotTitle=char(UniqueDirectionLabels(d));
    %             RowLabels=num2str(UniqueFreqs);
    %             PlotOLPupilSparklines(squeeze(AvgTimeSeries(:,d,:)),squeeze(AvgModelFits(:,d,:)),squeeze(SEMTimeSeries(:,d,:)),ErrorMultiplier,yRange, figSparkLinePlots(SubjectID), RowLabels, PlotTitle, d, length(UniqueDirections));
    %         end
    %     end % if statement catching edge case of one crossing of freq and direction
    %
    %
    %
    %     %% Plot avg pupil size by direction and frequency
    %
    %     PlotData=AvgMeanPupilByFxD;
    %     PlotErrors=StdMeanPupilByFxD*2;
    %
    %     PlotTitle=['Subject ' char(SubjectAnon) '/Mean Pupil Diameter [±2SD]'];
    %     PlotOLPupilTTFD(UniqueFreqs,UniqueDirectionLabels,PlotData,PlotErrors,PlotTitle,'Pupil diameter [mm]',[0 max(max(PlotData))+max(max(PlotErrors))],figMeanPupilByDirection,SubjectID,maxSubplots);
    %
    %
    %
    %     %     % Save out the data
    %     %
    %     %     if (SaveDataFlag==1)
    %     %
    %     %         currDir=pwd;
    %     %
    %     %         cd(ResultsFullPath)
    %     %         if ~exist(ResultsDirName)
    %     %             mkdir(ResultsDirName)
    %     %         end
    %     %
    %     %         cd(ResultsDirName);
    %     %
    %     %         dataFile = sprintf('%s-results.csv', Subjects{SubjectID});
    %     %         c = CSVFile(dataFile, true);
    %     %
    %     %         phase = AvgTimeSeriesPhases;
    %     %         amplitude = AvgTimeSeriesAmplitudes;
    %     %         phaseError = StdPhaseSplitHalfAvgTimeSeries;
    %     %         amplitudeError = StdAmplitudeSplitHalfAvgTimeSeries;
    %     %
    %     %         allDirections = repmat(UniqueDirectionLabels, length(UniqueFreqs), 1);
    %     %         c = c.addColumn('Direction', 's');
    %     %         c = c.setColumnData('Direction', allDirections(:));
    %     %
    %     %         c = c.addColumn('Frequency [Hz]', 'g');
    %     %         c = c.setColumnData('Frequency [Hz]', repmat(UniqueFreqs', length(UniqueDirectionLabels), 1));
    %     %
    %     %         c = c.addColumn('Pupil size mean [mm]', 'g');
    %     %         c = c.setColumnData('Pupil size mean [mm]', AvgMeanPupilByFxD(:));
    %     %
    %     %         c = c.addColumn('Pupil size error [1SD]', 'g');
    %     %         c = c.setColumnData('Pupil size error [1SD]', StdMeanPupilByFxD(:));
    %     %
    %     %         c = c.addColumn('Amplitude mean [proportion change]', 'g');
    %     %         c = c.setColumnData('Amplitude mean [proportion change]', amplitude(:));
    %     %
    %     %         c = c.addColumn('Phase mean (±pi)', 'g');
    %     %         c = c.setColumnData('Phase mean (±pi)', phase(:));
    %     %
    %     %         c = c.addColumn('Amplitude error [1SD]', 'g');
    %     %         c = c.setColumnData('Amplitude error [1SD]', amplitudeError(:));
    %     %
    %     %         c = c.addColumn('Phase error [1SD]', 'g');
    %     %         c = c.setColumnData('Phase error [1SD]', phaseError(:));
    %     %
    %     %         % Add contrasts (if config filenames were passed)
    %     %         if ~(isempty(configFileNames))
    %     %             for r = 1:length(receptorNames)
    %     %                 % Target contrasts
    %     %                 c = c.addColumn([receptorNames{r} ' contrast (target)'], 'g');
    %     %                 tmpData = (round(repmat(contrastTarget(:, r), 1, length(UniqueFreqs))*1000)/1000)';
    %     %                 c = c.setColumnData([receptorNames{r} ' contrast (target)'], tmpData(:));
    %     %
    %     %                 % Predicted contasts (with dark)
    %     %                 c = c.addColumn([receptorNames{r} ' contrast (predicted)'], 'g');
    %     %                 tmpData = (round(repmat(contrastPredicted(:, r), 1, length(UniqueFreqs))*1000)/1000)';
    %     %                 c = c.setColumnData([receptorNames{r} ' contrast (predicted)'], tmpData(:));
    %     %             end
    %     %         end
    %     %
    %     %         c.write;
    %     %
    %     %         cd(currDir);
    %     %     end % if SaveDataFlag
    %     %
    %     % end % subject loop
    %     %
    %     % fprintf('\nDone all subjects.\n'); % Notify user
    %     %
    %     %% Save out the plots
    %
    %     if (SavePlotsFlag==1)
    %
    %         fprintf('\nNow saving plots to disk.\n'); % Notify user
    %
    %         currDir=pwd;
    %
    %         cd(ResultsFullPath)
    %         if ~exist(ResultsDirName)
    %             mkdir(ResultsDirName)
    %         end
    %
    %         cd(ResultsDirName);
    %
    %
    %         % Save Avg Pupil size figure
    %         set(figMeanPupilByDirection, 'PaperPosition', [0 0 8 5]); %Position plot at left hand corner with width 5 and height 5.
    %         set(figMeanPupilByDirection, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
    %         saveas(figMeanPupilByDirection, 'mean_pupil_size_subjects.pdf', 'pdf');
    %
    %         % Save polar plots and sparkline plots
    %         for SubjectID=1:length(Subjects)
    %
    %             set(figSparkLinePlots(SubjectID), 'PaperPosition', [0 0 8 5]); %Position plot at left hand corner with width 5 and height 5.
    %             set(figSparkLinePlots(SubjectID), 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
    %             SubjectAnon = ['S' num2str(SubjectID)];
    %             saveas(figSparkLinePlots(SubjectID), ['sparkline_plot_subject_' char(SubjectAnon) '.pdf'], 'pdf');
    %
    %         end  % loop across subjects
    %
    %
    %         cd(currDir);
    %     end % if SavePlotsFlag
    %
    %
    %     pause off;
    
end % main
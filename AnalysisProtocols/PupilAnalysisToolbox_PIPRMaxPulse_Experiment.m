[~, userID] = system('whoami');
userID = strtrim(userID);


%% List of all the subjects and dates

Subjects = {'TEST_0004'};
Dates = {'072616'};

%% Params common to all components of the experiment

params = PupilAnalysisToolbox_GetDefaultParams;
params.BadNaNThreshold = 0.5;
params.PulseOnsetSecs = 1;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 1-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false;
params.full_trial_length = 17;
params.final_trial_length = 14;
params.minimum_trial_length = 14;
params.xLim = 14;

% Stimulus labels

oldLabels = {'PIPRMaxPulse-PulsePIPRBlue_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-PulsePIPRred_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-PulsePIPRRed_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundPIPR_45sSegment',...
             'PIPRMaxPulse-PulseMaxLMS_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundLMS_45sSegment',...
             'PIPRMaxPulse-PulseMaxMel_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundMel_45sSegment'};
newLabels = {'PIPRBlue', 'PIPRRed', 'PIPRRed', 'Background', 'MaxLMS', 'Background', 'MaxMel', 'Background'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulsePIPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Protocol={'PIPRMaxPulse_PulsePIPR'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulseLMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Protocol={'PIPRMaxPulse_PulseLMS'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulseLMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Protocol={'PIPRMaxPulse_PulseMel'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

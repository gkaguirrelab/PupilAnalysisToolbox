[~, userID] = system('whoami');
userID = strtrim(userID);


%% List of all the subjects and dates

Subjects = {'MELA_0039' 'MELA_0074'};
Dates = {'092716' '092716'};

%% Params common to all components of the experiment

params = PupilAnalysisToolbox_GetDefaultParams;
params.BadNaNThreshold = 0.5;
params.PulseOnsetSecs = 1;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 1-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = true;
params.full_trial_length = 17;
params.final_trial_length = 14;
params.minimum_trial_length = 14;
params.xLim = 14;
params.yLim = 0.7;

% Stimulus labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulsePIPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldLabels = {'PIPRMaxPulse-PulsePIPRBlue_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-PulsePIPRred_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundPIPR_45sSegment'};
newLabels = {'PIPRBlue', 'PIPRRed', 'Background'};

Protocol={'PIPRMaxPulse_PulsePIPR'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulseLMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldLabels = {'PIPRMaxPulse-PulseMaxLMS_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundLMS_45sSegment'};
newLabels = {'MaxLMS', 'Background'};

Protocol={'PIPRMaxPulse_PulseLMS'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulseLMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldLabels = {'PIPRMaxPulse-PulseMaxMel_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-BackgroundMel_45sSegment'};
newLabels = {'MaxMel', 'Background'};

Protocol={'PIPRMaxPulse_PulseMel'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/' Protocol{1}];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/' Protocol{1}];
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocol, newLabels, oldLabels, basePath, resultsPath);

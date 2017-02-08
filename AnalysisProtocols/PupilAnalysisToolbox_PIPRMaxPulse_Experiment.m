[~, userID] = system('whoami');
userID = strtrim(userID);


%% List of all the subjects and dates
% Master list
Subjects = {'MELA_0003' 'MELA_0026' 'MELA_0037' 'MELA_0038' 'MELA_0038' 'MELA_0039' 'MELA_0043' 'MELA_0049' 'MELA_0050' 'MELA_0071' 'MELA_0073' 'MELA_0074' 'MELA_0075' 'MELA_0077' 'MELA_0077' 'MELA_0078' 'MELA_0079' 'MELA_0080' 'MELA_0081' 'MELA_0082' 'MELA_0084' 'MELA_0085' 'MELA_0087' 'MELA_0088' 'MELA_0089' 'MELA_0090' 'MELA_0093' 'MELA_0093' 'MELA_0094' 'MELA_0096' 'MELA_0096' 'MELA_0097' 'MELA_0098' 'MELA_0100'};
Dates = {'110316' '112116' '120616' '020217' '092916' '092716' '092816' '102716' '110116' '100616' '102016' '092716' '111116' '020217' '110216' '100416' '121616' '102416' '101116' '113016' '102716' '110916' '112916' '111716' '111616' '112216' '020117' '120116' '121516' '020117' '120916' '010917' '121416' '121616'};
%%% DO NOT CHANGE THIS LIST APART FROM ADDING THE MOST RECENT SUBJECTS

%% Change the below
Subjects = { };
Dates = { };
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
params.yLim = 0.7;

% Stimulus labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulsePIPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldLabels = {'PIPRMaxPulse-PulsePIPRBlue_3s_MaxContrast17sSegment',...
             'PIPRMaxPulse-PulsePIPRRed_3s_MaxContrast17sSegment',...
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

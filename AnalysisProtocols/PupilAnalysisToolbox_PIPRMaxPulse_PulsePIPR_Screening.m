[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPRMaxPulse_PulsePIPR_Screening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Protocols={'PIPRMaxPulse_PulsePIPR_Screening'};
Subjects = {'TEST_0003'};
Dates = {'071316'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/PIPRMaxPulse_PulsePIPR_Screening'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/PIPRMaxPulse_PulsePIPR_Screening'];
newLabels = {'PIPRBlue', 'PIPRRed', 'PIPRRed'};
oldLabels = {'PIPRMaxPulse-PulsePIPRBlue_3s_MaxContrast17sSegment', 'PIPRMaxPulse-PulsePIPRRed_3s_MaxContrast17sSegment', 'PIPRMaxPulse-PulsePIPRred_3s_MaxContrast17sSegment'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocols, newLabels, oldLabels, basePath, resultsPath);


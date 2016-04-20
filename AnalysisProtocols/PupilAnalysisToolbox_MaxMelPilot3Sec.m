[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MaxMelLMSPilot_LMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.5;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds

Subjects = {'HERO_gka1'};
Protocols={'MaxMelLMSPilot_LMS'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MaxMelLMSPilot_LMS'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMelLMSPilot'];

oldLabels = {'MaxLMS400Pct-45sBackground', 'MaxLMS400Pct-45sPositivePulse3s'};
newLabels = {'LMSBackground', 'LMS400Pct'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params);

%%
Protocols={'MaxMelLMSPilot_Mel'};
params.PulseOnsetSecs = 4.5;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MaxMelLMSPilot_Mel'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMelLMSPilot'];

oldLabels = {'MaxMel400Pct-45sBackground', 'MaxMel400Pct-45sPositivePulse3s'};
newLabels = {'MelBackground', 'Mel400Pct'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params);
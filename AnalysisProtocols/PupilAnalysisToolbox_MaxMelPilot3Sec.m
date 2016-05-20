[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MaxMelLMSPilot_LMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 5;
params.PulseDurationSecs = 3;
params.yLim = 0.6;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds

Subjects = {'HERO_JAR_001'};
Protocols={'MelanopsinMRPupil_MaxMel400Pct3sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinMRPupil_MaxMel400Pct3sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMelLMSPilot'];

oldLabels = {'MaxMel400Pct-45sBackground', 'MaxMel400Pct-45sPositivePulse3s'};
newLabels = {'MelBackground', 'Mel400Pct'};
Data1 = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
%PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, [], []);

%%
Protocols={'MelanopsinMRPupil_MaxLMS400Pct3sPulse'};
params.PulseOnsetSecs = 5;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinMRPupil_MaxLMS400Pct3sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMelLMSPilot'];

oldLabels = {'MaxLMS400Pct-45sBackground', 'MaxLMS400Pct-45sPositivePulse3s'};
newLabels = {'LMSBackground', 'LMS400Pct'};
Data2 = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
%PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, [], []);


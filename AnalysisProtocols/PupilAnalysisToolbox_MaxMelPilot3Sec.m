[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_PIPR5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 5;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false;

Protocols={'MelanopsinMRPupil_MaxLMS400Pct3sPulse'};
Subjects = {'HERO_mxs1'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinMRPupil_MaxLMS400Pct3sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinMRPupil_MaxLMS400Pct3sPulse'];


newLabels = {'Background' 'LMS400Pct'}
oldLabels = {'MaxLMS400Pct-45sBackground', 'MaxLMS400Pct-45sPositivePulse3s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {[2 3]}, {[1 -1]});
%PupilAnalysisToolbox_SummarizeDataQuality(Subjects, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_PIPR5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 5;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false;

Protocols={'MelanopsinMRPupil_MaxMel400Pct3sPulse'};
Subjects = {'HERO_mxs1'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinMRPupil_MaxMel400Pct3sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinMRPupil_MaxMel400Pct3sPulse'];


newLabels = {'Background' 'Mel400Pct'}
oldLabels = {'MaxMel400Pct-45sBackground', 'MaxMel400Pct-45sPositivePulse3s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {[2 3]}, {[1 -1]});
%PupilAnalysisToolbox_SummarizeDataQuality(Subjects, resultsPath);
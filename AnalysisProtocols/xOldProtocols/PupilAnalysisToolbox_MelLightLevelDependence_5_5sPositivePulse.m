[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sPulse_ND00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G092815AxND00' 'J092915RxND00' 'M092515SxND00'};
Protocols={'MelLightLevelDependence5sPulse_ND00'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sPulse_ND00'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sPulse_ND00'];

newLabels = {'Background', 'LMS+', 'Mel+',};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sPulse_ND05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G100515AxND05' 'J100615RxND05' 'M100615SxND05'};
Protocols={'MelLightLevelDependence5sPulse_ND05'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sPulse_ND05'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sPulse_ND05'];

newLabels = {'Background', 'LMS+', 'Mel+',};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sPulse_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G092815AxND10' 'J092915RxND10' 'M092515SxND10'};
Protocols={'MelLightLevelDependence5sPulse_ND10'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sPulse_ND10'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sPulse_ND10'];


newLabels = {'Background', 'LMS+', 'Mel+',};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sPulse_ND15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G100515AxND15' 'J100615RxND15' 'M100615SxND15'};
Protocols={'MelLightLevelDependence5sPulse_ND15'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sPulse_ND15'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sPulse_ND15'];

newLabels = {'Background', 'LMS+', 'Mel+',};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sPulse_ND20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G092815AxND20' 'J092915RxND20' 'M092515SxND20'};
Protocols={'MelLightLevelDependence5sPulse_ND20'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sPulse_ND20'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sPulse_ND20'];

newLabels = {'Background', 'LMS+', 'Mel+',};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sBipolarPulse_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G100815AxND10' 'J100715RxND10' 'M100915SxND10'};
Protocols={'MelLightLevelDependence5sBipolarPulse_ND10'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sBipolarPulse_ND10'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sBipolarPulse_ND10'];

newLabels = {'Background', 'LMS+', 'LMS-', 'Mel+', 'Mel-'};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'LMSDirectedNulled-45sNegativePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sNegativePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sBipolarPulse_ND15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G100815AxND15' 'J100715RxND15' 'M100915SxND15'};
Protocols={'MelLightLevelDependence5sBipolarPulse_ND15'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sBipolarPulse_ND15'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sBipolarPulse_ND15'];

newLabels = {'Background', 'LMS+', 'LMS-', 'Mel+', 'Mel-'};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5s', 'LMSDirectedNulled-45sNegativePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5s', 'MelanopsinDirectedPenumbralIgnoreNulled-45sNegativePulse5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sBipolarPulseConeNoise_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

Subjects = {'G100815AxND10' 'M100915SxND10' 'J100715RxND10'};
Protocols={'MelLightLevelDependence5sBipolarPulseConeNoise_ND10'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/xCompleted/CommonBackgroundPupil_ARVO2016/MelLightLevelDependence5sBipolarPulseConeNoise_ND10'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sBipolarPulseConeNoise_ND10'];

newLabels = {'Background', 'LMS+', 'LMS-', 'Mel+', 'Mel-'};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5sConeNoise', 'LMSDirectedNulled-45sNegativePulse5sConeNoise', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5sConeNoise', 'MelanopsinDirectedPenumbralIgnoreNulled-45sNegativePulse5sConeNoise'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {});
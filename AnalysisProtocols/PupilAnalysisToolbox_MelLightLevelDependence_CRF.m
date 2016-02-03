%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnoreNegativePulseConeNoiseCRF_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'J100715RxND10' 'G100815AxND10' 'M100915SxND10'};
Protocols={'MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnoreNegativePulseConeNoiseCRF_ND10'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnoreNegativePulseConeNoiseCRF_ND10';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sCRF_ND10';

newLabels = {'Mel_-09pct', 'Mel_-18pct', 'Mel_-36pct'};
oldLabels = {'MelanopsinDirectedPenumbralIgnoreNulled9Pct-45sNegativePulse5sConeNoiseCRF' 'MelanopsinDirectedPenumbralIgnoreNulled18Pct-45sNegativePulse5sConeNoiseCRF' 'MelanopsinDirectedPenumbralIgnoreNulled36Pct-45sNegativePulse5sConeNoiseCRF'};    
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnorePositivePulseConeNoiseCRF_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'J100715RxND10' 'G100815AxND10' 'M100915SxND10'};
Protocols={'MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnorePositivePulseConeNoiseCRF_ND10'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelLightLevelDependence5sMelanopsinDirectedPenumbralIgnorePositivePulseConeNoiseCRF_ND10';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sCRF_ND10';

newLabels = {'Mel_+09pct', 'Mel_+18pct', 'Mel_+36pct'};
oldLabels = {'MelanopsinDirectedPenumbralIgnoreNulled9Pct-45sPositivePulse5sConeNoiseCRF' 'MelanopsinDirectedPenumbralIgnoreNulled18Pct-45sPositivePulse5sConeNoiseCRF' 'MelanopsinDirectedPenumbralIgnoreNulled36Pct-45sPositivePulse5sConeNoiseCRF'};    
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sLMSDirectedNegativePulseConeNoiseCRF_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'J100715RxND10' 'G100815AxND10' 'M100915SxND10'};
Protocols={'MelLightLevelDependence5sLMSDirectedNegativePulseConeNoiseCRF_ND10'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelLightLevelDependence5sLMSDirectedNegativePulseConeNoiseCRF_ND10';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sCRF_ND10';

newLabels = {'LMS_-09pct', 'LMS_-18pct', 'LMS_-36pct'};
oldLabels = {'LMSDirectedNulled9Pct-45sNegativePulse5sConeNoiseCRF' 'LMSDirectedNulled18Pct-45sNegativePulse5sConeNoiseCRF' 'LMSDirectedNulled36Pct-45sNegativePulse5sConeNoiseCRF'};    
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelLightLevelDependence5sLMSDirectedPositivePulseConeNoiseCRF_ND10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'J100715RxND10' 'G100815AxND10' 'M100915SxND10'};
Protocols={'MelLightLevelDependence5sLMSDirectedPositivePulseConeNoiseCRF_ND10'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelLightLevelDependence5sLMSDirectedPositivePulseConeNoiseCRF_ND10';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelLightLevelDependence5sCRF_ND10';

newLabels = {'LMS_+09pct', 'LMS_+18pct', 'LMS_+36pct'};
oldLabels = {'LMSDirectedNulled9Pct-45sPositivePulse5sConeNoiseCRF' 'LMSDirectedNulled18Pct-45sPositivePulse5sConeNoiseCRF' 'LMSDirectedNulled36Pct-45sPositivePulse5sConeNoiseCRF'};    
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);


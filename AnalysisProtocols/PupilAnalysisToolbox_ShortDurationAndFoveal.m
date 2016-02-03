[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 1;

Subjects = {'M012216S' 'J012216R' 'G012216A'};
Protocols={'MelanopsinStepsFovealControlShortDuration1_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinStepsFovealControlShortDuration1_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinStepsFovealControlShortDuration1_5sPulse'];

newLabels = {'Background', 'LMS+', 'Mel+', 'Cone noise'};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse1_5sConeNoise', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse1_5sConeNoise', 'ConeNoiseOnly-45s'};
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'M012216S' 'J012216R' 'G012216A'};
Protocols={'MelanopsinStepsFovealControlShortDuration5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinStepsFovealControlShortDuration5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinStepsFovealControlShortDuration5_5sPulse'];

newLabels = {'Background', 'LMS+', 'Mel+', 'Cone noise'};
oldLabels = {'Background-60s', 'LMSDirectedNulled-45sPositivePulse5_5sConeNoise', 'MelanopsinDirectedPenumbralIgnoreNulled-45sPositivePulse5_5sConeNoise', 'ConeNoiseOnly-45s'};
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mel5_5sPulseFoveal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;

Subjects = {'J012216R_foveal'};
Protocols={'Mel5_5sPulseFoveal'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/Mel5_5sPulseFoveal'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/Mel5_5sPulseFoveal'];
newLabels = {'Background', 'LMS+', 'Mel+', 'Cone noise'};
oldLabels = {'Background-60s', 'LMSDirectedNulledFoveal-45sPositivePulse5_5sConeNoise', 'MelanopsinDirectedPenumbralIgnoreNulledFoveal-45sPositivePulse5_5sConeNoise', 'ConeNoiseOnly-45s'};
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);

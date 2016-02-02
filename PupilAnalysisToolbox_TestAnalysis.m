%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams
params.StepDurSecs = 1;

Subjects = {'M012216S' 'J012216R' 'G012216A'};
Protocols={'MelanopsinStepsFovealControlShortDuration1_5sPulse'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinStepsFovealControlShortDuration1_5sPulse';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinStepsFovealControlShortDuration1_5sPulse';

newLabels = {'Background', 'LMS+', 'Mel+', 'Cone noise'};
oldLabels = {'Background', 'LMSDirectedNulled', 'MelanopsinDirectedPenumbralIgnoreNulled', 'ConeNoiseOnly'};
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams
params.StepDurSecs = 5;

Subjects = {'M012216S' 'J012216R' 'G012216A'};
Protocols={'MelanopsinStepsFovealControlShortDuration5_5sPulse'};

basePath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinStepsFovealControlShortDuration5_5sPulse';
resultsPath = '/Users/pupillab/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinStepsFovealControlShortDuration5_5sPulse';

newLabels = {'Background', 'LMS+', 'Mel+', 'Cone noise'};
oldLabels = {'Background', 'LMSDirectedNulled', 'MelanopsinDirectedPenumbralIgnoreNulled', 'ConeNoiseOnly'};
PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
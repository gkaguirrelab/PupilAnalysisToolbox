[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;
params.yLim = 0.8;

Subjects = {'TEST_0001'};
Protocols={'SilentSubstitutionPIPR_PIPR5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_PIPR5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionPIPR_PIPR5_5sPulse'];

newLabels = {'Background-60s', 'Background-45s', 'PIPRBlue', 'PIPRRed'};
oldLabels = {'SilentSubstitutionPIPRBackgroundPIPR-60s', 'SilentSubstitutionPIPRBackgroundPIPR-45s', 'SilentSubstitutionPIPRBlue-45sPositivePulse5_5s', 'SilentSubstitutionPIPRRed-45sPositivePulse5_5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.yLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;
params.yLim = 0.4;

Subjects = {'TEST_0001'};
Protocols={'SilentSubstitutionPIPR_SS5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_SS5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionSS_PIPR5_5sPulse'];

newLabels = {'Background-60s', 'Cone noise', 'LMS+', 'Mel+'};
oldLabels = {'SilentSubstitutionPIPRBackgroundSS-60s', 'SilentSubstitutionPIPRConeNoiseOnly-45s', 'SilentSubstitutionPIPRLMSDirected-45sPositivePulse5_5sConeNoise', 'SilentSubstitutionPIPRMelanopsinDirected-45sPositivePulse5_5sConeNoise'};

Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.yLim);

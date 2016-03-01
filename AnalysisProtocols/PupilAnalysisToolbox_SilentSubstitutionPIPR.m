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

basePath = ['~/Desktop/SilentSubstitutionPIPR_PIPR5_5sPulse'];
resultsPath = ['~/Desktop/SilentSubstitutionPIPR_PIPR5_5sPulse'];

newLabels = {'Background-60s', 'Background-45s', 'PIPRBlue', 'PIPRRed'};
oldLabels = {'SilentSubstitutionPIPRBackgroundPIPR-60s', 'SilentSubstitutionPIPRBackgroundPIPR-45s', 'SilentSubstitutionPIPRBlue-45sPositivePulse5_5s', 'SilentSubstitutionPIPRRed-45sPositivePulse5_5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.yLim);

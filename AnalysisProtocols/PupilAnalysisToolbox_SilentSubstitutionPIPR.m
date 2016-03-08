[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_PIPR5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;
params.yLim = 0.8;

Subjects = {'MELA_0043'};
Protocols={'SilentSubstitutionPIPR_PIPR5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_PIPR5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionPIPR_PIPR5_5sPulse'];

newLabels = {'Background', 'PIPRBlue', 'PIPRRed'};
oldLabels = {'SilentSubstitutionPIPRBackgroundPIPR-45s', 'SilentSubstitutionPIPRBlue-45sPositivePulse5_5s', 'SilentSubstitutionPIPRRed-45sPositivePulse5_5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.xLim, params.yLim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_SS5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;
params.yLim = 0.8;

Subjects = {'MELA_0043'};
Protocols={'SilentSubstitutionPIPR_SS5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_SS5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionSS_PIPR5_5sPulse'];

newLabels = {'ConeNoise', 'LMS+', 'Mel+'};
oldLabels = {'SilentSubstitutionPIPRBackgroundSS-45s', 'SilentSubstitutionPIPRLMSDirected-45sPositivePulse5_5sConeNoise', 'SilentSubstitutionPIPRMelanopsinDirected-45sPositivePulse5_5sConeNoise'};

Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.xLim, params.yLim);

[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 1.5;
params.yLim = 0.8;

Subjects = {'TEST_0002'};
Protocols={'MaxMel_1_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MaxMel_1_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMel_1_5sPulse'];

oldLabels = {'MaxMelBackground-45s', 'MaxMelPos-45sPositivePulse1_5s'};
newLabels = {'ConeNoise' 'Mel+'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.yLim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 1.5;
params.yLim = 0.8;

Subjects = {'TEST_0001'};
Protocols={'MaxMel_MaxMelIntermixedPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MaxMel_MaxMelIntermixedPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMel_MaxMelIntermixedPulse'];

oldLabels = {'MaxMelBackground-45s', 'MaxMelPos-45sPositivePulse5_5s', 'MaxMelPos-45sPositivePulse1_5s'};
newLabels = {'ConeNoise' 'Mel+ 5.5 sec' 'Mel+ 1.5 sec'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.yLim);
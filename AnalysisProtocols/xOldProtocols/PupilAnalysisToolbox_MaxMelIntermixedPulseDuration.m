[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinStepsFovealControlShortDuration1_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.StepDurSecs = 5.5;
params.yLim = 0.8;

Subjects = {'HERO_gka1'}; % need to fix for mxs1
Protocols={'MaxMel_MaxMelIntermixedPulseDuration'};


basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MaxMel_MaxMelIntermixedPulseDuration'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MaxMel_MaxMelIntermixedPulseDuration'];

oldLabels = {'MaxMelConeNoise-45s' 'MaxMelPos-45sPositivePulse1_5s' 'MaxMelPos-45sPositivePulse2_5s' 'MaxMelPos-45sPositivePulse3_5s' 'MaxMelPos-45sPositivePulse4_5s' 'MaxMelPos-45sPositivePulse5_5s'};
newLabels = {'ConeNoise' 'MaxMel_1_5s' 'MaxMel_2_5s' 'MaxMel_3_5s' 'MaxMel_4_5s' 'MaxMel_5_5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params.xLim, params.yLim);

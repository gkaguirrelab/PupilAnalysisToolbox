[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinPupilMaxLMS_3sPulse17sSegment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.BadNaNThreshold = 0.5;
params.PulseOnsetSecs = 1;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 1-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false;
params.full_trial_length = 17;
params.final_trial_length = 14;
params.minimum_trial_length = 14;
params.xLim = 14;
Protocols={'MelanopsinPupilMaxLMS_3sPulse17sSegment'};
Subjects = {'HERO_asb1' 'HERO_mxs1' 'HERO_gka1'};
Dates = {'061516' '061616' '061716'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinPupilMaxLMS_3sPulse17sSegment'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinPupilMaxLMS_3sPulse17sSegment'];
newLabels = {'LMS400Pct'}
oldLabels = {'MelanopsinPupilMaxLMS-PulseMaxLMS_3s_MaxContrast17sSegment'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocols, newLabels, oldLabels, basePath, resultsPath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MelanopsinPupilMaxLMS_3sPulse17sSegment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.BadNaNThreshold = 0.5;
params.PulseOnsetSecs = 1;
params.PulseDurationSecs = 3;
params.meanCenterWindow = [0 1-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false;
params.full_trial_length = 17;
params.final_trial_length = 14;
params.minimum_trial_length = 14;
params.xLim = 14;
Protocols={'MelanopsinPupilMaxMel_3sPulse17sSegment'};
Subjects = {'HERO_asb1' 'HERO_mxs1' 'HERO_gka1'};
Dates = {'061516' '061616' '061716'};
basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/MelanopsinPupilMaxMel_3sPulse17sSegment'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinPupilMaxMel_3sPulse17sSegment'];
newLabels = {'Mel400Pct'}
oldLabels = {'MelanopsinPupilMaxMel-PulseMaxMel_3s_MaxContrast17sSegment'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Dates, Protocols, newLabels, oldLabels, basePath, resultsPath);
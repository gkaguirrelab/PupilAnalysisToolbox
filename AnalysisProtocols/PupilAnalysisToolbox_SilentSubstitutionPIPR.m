[~, userID] = system('whoami');
userID = strtrim(userID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_PIPR5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

% This is a list of all the subjects

% Subjects = {'MELA_0001' 'MELA_0002' 'MELA_0017' 'MELA_0026' 'MELA_0028' ...
%    'MELA_0038' 'MELA_0039' 'MELA_0043' 'MELA_0044' 'MELA_0045' ...
%    'MELA_0046' 'MELA_0047' 'MELA_0049' 'MELA_0050' 'MELA_0051' ...
%    'MELA_0052' 'MELA_0053' 'MELA_0054' 'MELA_0055' 'MELA_0057' ...
%    'MELA_0058' 'MELA_0061' 'MELA_0062' 'MELA_0063' 'MELA_0065' ...
%    'MELA_0067' 'MELA_0068' 'MELA_0069' 'MELA_0070' 'MELA_0071' }

% This is a list of only the "bad" subjects

Subjects = {'MELA_0017' 'MELA_0026' 'MELA_0028' ...
   'MELA_0039' 'MELA_0044'  ...
   'MELA_0046' 'MELA_0047' 'MELA_0051' ...
   'MELA_0052' 'MELA_0054' 'MELA_0055' ...
   'MELA_0058' 'MELA_0063'  ...
   'MELA_0067' 'MELA_0070'  }

% This is a list of only the "good" subjects (who have had their data
%   merged and passed data quality check)

% Subjects = {'MELA_0001' 'MELA_0002' ...
%    'MELA_0038' 'MELA_0043' 'MELA_0045' ...
%    'MELA_0049' 'MELA_0050'  ...
%    'MELA_0053' 'MELA_0057' ...
%    'MELA_0061' 'MELA_0062' 'MELA_0065' ...
%    'MELA_0068' 'MELA_0069' 'MELA_0071' };

% This is a subject listing that could be used to examine the data from a
% single subject and check data quality

% Subjects = {'MELA_0001' 'MELA_0002'}

Protocols={'SilentSubstitutionPIPR_PIPR5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_PIPR5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionPIPR_PIPR5_5sPulse'];

newLabels = {'Background', 'PIPRBlue', 'PIPRRed'};
oldLabels = {'SilentSubstitutionPIPRBackgroundPIPR-45s', 'SilentSubstitutionPIPRBlue-45sPositivePulse5_5s', 'SilentSubstitutionPIPRRed-45sPositivePulse5_5s'};
Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {[2 3]}, {[1 -1]});
PupilAnalysisToolbox_SummarizeDataQuality(Subjects, resultsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SilentSubstitutionPIPR_SS5_5sPulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = PupilAnalysisToolbox_GetDefaultParams;
params.PulseOnsetSecs = 4.75;
params.PulseDurationSecs = 5.5;
params.meanCenterWindow = [0 params.PulseOnsetSecs-1/params.sampling_frequency]; % In seconds
params.TrialInspectorFlag = false; % trial-by-trial turned off

% This is a list of all the subjects

% Subjects = {'MELA_0001' 'MELA_0002' 'MELA_0017' 'MELA_0026' 'MELA_0028' ...
%    'MELA_0038' 'MELA_0039' 'MELA_0043' 'MELA_0044' 'MELA_0045' ...
%    'MELA_0046' 'MELA_0047' 'MELA_0049' 'MELA_0050' 'MELA_0051' ...
%    'MELA_0052' 'MELA_0053' 'MELA_0054' 'MELA_0055' 'MELA_0057' ...
%    'MELA_0058' 'MELA_0061' 'MELA_0062' 'MELA_0063' 'MELA_0065' ...
%    'MELA_0067' 'MELA_0068' 'MELA_0069' 'MELA_0071' }

% This is a list of only the "bad" subjects

Subjects = {'MELA_0017' 'MELA_0026' 'MELA_0028' ...
   'MELA_0039' 'MELA_0044'  ...
   'MELA_0046' 'MELA_0047' 'MELA_0051' ...
   'MELA_0052' 'MELA_0054' 'MELA_0055' ...
   'MELA_0058' 'MELA_0063'  ...
   'MELA_0067'  }

% This is a list of only the "good" subjects (who have had their data
%   merged and passed data quality check)

% Subjects = {'MELA_0001' 'MELA_0002' ...
%    'MELA_0038' 'MELA_0043' 'MELA_0045' ...
%    'MELA_0049' 'MELA_0050'  ...
%    'MELA_0053' 'MELA_0057' ...
%    'MELA_0061' 'MELA_0062' 'MELA_0065' ...
%    'MELA_0068' 'MELA_0069' 'MELA_0071' };

% This is a subject listing that could be used to examine the data from a
% single subject and check data quality

% Subjects = {'MELA_0071' }

Protocols={'SilentSubstitutionPIPR_SS5_5sPulse'};

basePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/SilentSubstitutionPIPR_SS5_5sPulse'];
resultsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/SilentSubstitutionSS_PIPR5_5sPulse'];

newLabels = {'ConeNoise', 'LMS+', 'Mel+'};
oldLabels = {'SilentSubstitutionPIPRBackgroundSS-45s', 'SilentSubstitutionPIPRLMSDirected-45sPositivePulse5_5sConeNoise', 'SilentSubstitutionPIPRMelanopsinDirected-45sPositivePulse5_5sConeNoise'};

Data = PupilAnalysisToolbox_PulseSequentialTrialAnalysis(params, Subjects, Protocols, newLabels, oldLabels, basePath, resultsPath);
PupilAnalysisToolbox_PlotTimeSeriesAndMean(Data, Subjects, resultsPath, params, {[3 2]}, {[1 -1]});
PupilAnalysisToolbox_SummarizeDataQuality(Subjects, resultsPath);

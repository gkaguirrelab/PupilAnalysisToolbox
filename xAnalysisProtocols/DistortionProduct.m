
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANTS AND CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

StimOnsetDelay=0; % This is the delay in msecs between the onset of
% the eye tracking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.


TrialInspectorFlag=0; % If set, plot each time series

% Generally, only one or the other of the two harmonic flags should be set

HarmonicModelFlag=1; % If set, model in the sparkline plots the harmonic
HarmonicTestFlag=0; % If set, generate all response plots assuming harmonic stimulation

adapt_length = 60; % Duration in seconds of the initial background adaptation peiod

sampling_frequency = 20;  % The samples per second of our final data vectors
full_trial_length = 45;  % The length in seconds of the full trial
%   acquisition
final_trial_length = 40; % The length of each data trial after processing,
%   taken from the tail of each full data trial.
%    That is, we discard some data from the start of each trial.
minimum_length_trial=42; % Any trial shorter than this will be discarded

StutterErrorFlag = 0;     % If these are data from the initial TTF4D
% experiment, use a modified version of the LeastSquaresSpectralFit routine
% that accounts for the tiny error in the stimulus sinusoidal modulation
spike_remover_params(1)=10; % The window size to use for spike removal
spike_remover_params(2)=0.2; % The max acceptable proportion change in the window
spike_remover_params(3)=2; % The max acceptable change in SD units in the window
sgolay_span = 20;         % Parameters of the SGolay interpolation prior
sgolay_polynomial = 7;    %   to interpolation
BadPercentChangeThreshold = .40;       % Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is
% certainly noise
BadNanThreshold=0.15; % If more than this proportion
% of a trial time series is composed of NaNs, discard the trial

% The setting of the Protocols variable will determine which
%  sessions will be concatenated together in the analysis

Subjects = {'s004Brimonidine'};

RelabelSubjectsFlag=0;

Protocols={'PupillometryDistortionProduct'};

configFileNames=[];
% NOTE: In later experiments other than TTF4D, this has to be set to 51.
whichSttingIndexToValidate = 51;




% If set to one, create a synthetic direction in the results which is the
% sum of the LM and Mel mods. Note that later this flag leads to hard-coded
% behavior in which the 2nd and 3rd modulation directions are assumed to be
% the Melanopsin and LM directions. This could break if other directions
% are included besides the standard four.

% The SequentialTriaAnalysis code detects the case of this flag set, as
% well as the OPN4ScaleFlag to ensure proper behavior

SyntheticIsoFlag = 0;

% When the OPN4ScaleFlag is set, scaled polar plots are created for each
% subject using the SyntheticIsoValues.

OPN4ScaleFlag = 0;

%% Handle plotting functions across subjects. This includes setting a look-
%% up table to replace the modulation directions with more succint labels

StackPlotFlag=1; % When set, multiple-subject results will be overlayed in
% a single plot

SaveDataFlag = 1;
SavePlotsFlag = 1;
RelabelDirections = 1;

%% Are we in v2 of the pupil stuff?
OLFlickerSensitivityFlag = 1;
analyzeGazeData = false;

% Add the OLPupilDiameter folder to the path. Different conventions.
addpath(genpath('/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter'));
addpath(genpath('/Users/Shared/Matlab/Experiments/OneLight/OLFlickerSensitivity'));


newLabels = {'Isochromatic'};

oldLabels = {'LMSDirected'};

ResultsDirName = 'PupillometryDistortionProduct';


TimeSeries=SequentialTrialAnalysis(...
    sampling_frequency,...
    full_trial_length,...
    final_trial_length,...
    adapt_length,...
    minimum_length_trial,...
    StutterErrorFlag,...
    spike_remover_params,...
    sgolay_span,...
    sgolay_polynomial,...
    BadPercentChangeThreshold,...
    BadNanThreshold,...
    Subjects,...
    Protocols,...
    SyntheticIsoFlag,...
    SaveDataFlag,...
    SavePlotsFlag,...
    RelabelDirections,...
    newLabels,...
    oldLabels,...
    ResultsDirName,...
    configFileNames,...
    whichSttingIndexToValidate,...
    TrialInspectorFlag,...
    HarmonicModelFlag,...
    HarmonicTestFlag,...
    StimOnsetDelay,...
    RelabelSubjectsFlag,...
    StackPlotFlag,...
    OPN4ScaleFlag,...
    OLFlickerSensitivityFlag,...
    analyzeGazeData);




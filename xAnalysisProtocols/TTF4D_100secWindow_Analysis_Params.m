


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANTS AND CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

OLFlickerSensitivityFlag=0; % Indicates that we are in the old-school
% sequential trial analysis state


StimOnsetDelay=1031; % This is the delay in msecs between the onset of 
% the eye trcking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.

TrialInspectorFlag=0; % If set, plot each time series

% Generally, only one or the other of the two harmonic flags should be set

HarmonicModelFlag=1; % If set, model in the sparkline plots the harmonic
HarmonicTestFlag=0; % If set, generate all response plots assuming harmonic stimulation

adapt_length = 300; % Duration in seconds of the initial background adaptation peiod

sampling_frequency = 20;  % The samples per second of our final data vectors
full_trial_length = 120;  % The length in seconds of the full trial
%   acquisition
final_trial_length = 100; % The length of each data trial after processing,
%   taken from the tail of each full data trial.
%    That is, we discard some data from the start of each trial.
minimum_length_trial=116; % Any trial shorter than this will be discarded

StutterErrorFlag = 1;     % If these are data from the initial TTF4D
% experiment, use a modified version of the LeastSquaresSpectralFit routine
% that accounts for the tiny error in the stimulus sinusoidal modulation
spike_remover_params(1)=10; % The window size to use for spike removal
spike_remover_params(2)=0.2; % The max acceptable proportion change in the window
spike_remover_params(3)=3; % The max acceptable change in SD units in the window

sgolay_span = 20;         % Parameters of the SGolay interpolation prior
sgolay_polynomial = 7;    %   to interpolation
BadPercentChangeThreshold = .40;       % Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is
% certainly noise
BadNanThreshold=0.10; % If more than this proportion
% of a trial time series is composed of NaNs, discard the trial


% The setting of the Protocols variable will determine which
%  sessions will be concatenated together in the analysis

Subjects={'Sandeep','Geoff'};

RelabelSubjectsFlag=1;

Protocols={'TTF4D_Iso_Peri_Dilate_5mmPinhole',...
    'TTF4D_Mel_Peri_Dilate_5mmPinhole',...
    'TTF4D_LM_Peri_Dilate_5mmPinhole',...
    'TTF4D_S_Peri_Dilate_5mmPinhole',...
    'TTF4D_Subset_Mel_Peri_Dilate_5mmPinhole',...
    'TTF4D_Subset_Iso_Peri_Dilate_5mmPinhole',...
    'TTF4D_Subset_Mel_Peri_Dilate_5mmPinhole',...
    'TTF4D_Subset_S_Peri_Dilate_5mmPinhole'};

% configFileNames = {'TTF4D-Isochromatic-OLEyeTrackerLongCable.cfg' ;
%     'TTF4D-LMDirected-OLEyeTrackerLongCable.cfg' ; ...
%     'TTF4D-MelanopsinDirected-OLEyeTrackerLongCable.cfg' ; ...
%     'TTF4D-SIsolating-OLEyeTrackerLongCable.cfg'};
configFileNames=[];

% NOTE: In later experiments other than TTF4D, this has to be set to 51.
whichSttingIndexToValidate = 50;

% If set to one, create a synthetic direction in the results which is the
% sum of the LM and Mel mods. Note that later this flag leads to hard-coded
% behavior in which the 2nd and 3rd modulation directions are assumed to be
% the Melanopsin and LM directions. This could break if other directions
% are included besides the standard four.

% The SequentialTriaAnalysis code detects the case of this flag set, as
% well as the OPN4ScaleFlag to ensure proper behavior

SyntheticIsoFlag = 1;

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

newLabels = {'Iso' ; 'LM' ; 'Mel' ; 'S'};
oldLabels = {'TTF4D-Isochromatic' ; 'TTF4D-LMDirected' ;...
    'TTF4D-MelanopsinDirected' ; 'TTF4D-SIsolating'};

ResultsDirName = 'TTF4D_100secWindow';


% Add the OLPupilDiameter folder to the path. Different conventions.
if isdir('/Users/Shared/Matlab/Experiments/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OLPupilDiameter';
elseif isdir('/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/');
    basePath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter';
end
addpath(genpath(basePath));


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
    OLFlickerSensitivityFlag);

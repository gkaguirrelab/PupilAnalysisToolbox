
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANTS AND CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

 % This is the delay in msecs between the onset of 
% the eye trcking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.
StimOnsetDelay=1031;

% Do we have different contrast levels?
contrastAnalysis = 1;

% If set, plot each time series
TrialInspectorFlag=0; 

% Number of sessions
num_session = 20;

% Set to 1 if we want to compute the MTF of the first harmonic
HarmonicModelFlag=0;

% Duration in seconds of the initial background adaptation peiod
adapt_length = 300; 

% The samples per second of our final data vectors
sampling_frequency = 50;  

% The length in seconds of the full trial acquisition
full_trial_length = 50; 

% The length of each data trial after processing,
% That is, we discard some data from the start and the end of each trial.
final_trial_length = 40; 

% Any trial shorter than this will be discarded
minimum_length_trial=40; 

 % If these are data from the initial TTF4D
% experiment, use a modified version of the LeastSquaresSpectralFit routine
% that accounts for the tiny error in the stimulus sinusoidal modulation
StutterErrorFlag = 0;    

% Parameters for spike removal algorithm
spike_remover_params(1)=10; % The window size
spike_remover_params(2)=0.2; % The max acceptable proportion change in the window
spike_remover_params(3)=2; % The max acceptable change in SD units in the window

% Parameters of Savitzky-Golay filter including the window size (sgolay_span)
% and the degree of smoothing polynomial (sgolay_polynomial)
sgolay_span = 20;         
sgolay_polynomial = 7;    

% Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is
% certainly noise
BadPercentChangeThreshold = .40;  

% If more than this proportion
% of a trial time series is composed of NaNs, discard the trial
BadNanThreshold=0.10; 

% The setting of the Protocols variable will determine which
%  sessions will be concatenated together in the analysis

% Subject list
Subjects = {'M112813S'};

RelabelSubjectsFlag=0;

Protocols={'ContrastEnvelopePupillometry'};

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
OPN4ScaleFlag = 1;

%% Handle plotting functions across subjects. This includes setting a look-
%% up table to replace the modulation directions with more succint labels

StackPlotFlag=1; % When set, multiple-subject results will be overlayed in
% a single plot

SaveDataFlag = 1;
SavePlotsFlag = 1;
RelabelDirections = 1;

newLabels = {'LM' ; 'Mel' ; 'S' ; 'Isochromatic' ; 'RodDirected'};
oldLabels = {'Modulation-LMDirected-50sContrastModulation-26' ; 'Modulation-MelanopsinDirected-50sContrastModulation-26' ;...
    'Modulation-SDirected-50sContrastModulation-26' ; 'Modulation-Isochromatic-50sContrastModulation-26' ; 'Modulation-RodDirectedHigh-50sContrastModulation-26'};


ResultsDirName = 'Envelope_0.5Hz40secWindow';


TimeSeries=SequentialTrialAnalysisEnvelope(...
    num_session,...
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
    contrastAnalysis,...
    RelabelDirections,...
    ResultsDirName,...
    configFileNames,...
    whichSttingIndexToValidate,...
    TrialInspectorFlag,...
    HarmonicModelFlag,...
    StimOnsetDelay,...
    RelabelSubjectsFlag);





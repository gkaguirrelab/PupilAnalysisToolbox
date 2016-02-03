


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANTS AND CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

StimOnsetDelay=1031; % This is the delay in msecs between the onset of 
% the eye trcking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.


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
SpikeFilterFlag = 1;      % Do we apply the spike filter routine?
sgolay_span = 20;         % Parameters of the SGolay interpolation prior
sgolay_polynomial = 7;    %   to interpolation
BadThreshold = .15;       % Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is
% certainly noise


% The setting of the Protocols variable will determine which
%  sessions will be concatenated together in the analysis

Subjects={'Test_Phase'};
Protocols={'TTF4D_Iso_Peri_Dilate_5mmPinhole'};

configFileNames = {'TTF4D-Isochromatic-OLEyeTrackerLongCable.cfg'};
% NOTE: In later experiments other than TTF4D, this has to be set to 51.
whichSttingIndexToValidate = 50;

% If set to one, create a synthetic direction in the results which is the
% sum of the LM and Mel mods. Note that later this flag leads to hard-coded
% behavior in which the 2nd and 3rd modulation directions are assumed to be
% the Melanopsin and LM directions. This could break if other directions
% are included besides the standard four.

SyntheticIsoFlag = 0;

% Handle plotting functions across subjects. This includes setting a look-
% up table to replace the modulation directions with more succint labels

SaveDataFlag = 0;
SavePlotsFlag = 0;
RelabelDirections = 0;

newLabels = {'Iso' ; 'LM' ; 'Mel' ; 'S'};
oldLabels = {'TTF4D-Isochromatic' ; 'TTF4D-LMDirected' ;...
    'TTF4D-MelanopsinDirected' ; 'TTF4D-SIsolating'};

ResultsDirName = 'RandPhaseTest_100secWindow';


TimeSeries=SequentialTrialAnalysis(...
    sampling_frequency,...
    full_trial_length,...
    final_trial_length,...
    minimum_length_trial,...
    StutterErrorFlag,...
    SpikeFilterFlag,...
    sgolay_span,...
    sgolay_polynomial,...
    BadThreshold,...
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
    whichSttingIndexToValidate);


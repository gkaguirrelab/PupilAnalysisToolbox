function params = PupilAnalysisToolbox_GetDefaultParams

params.StimOnsetDelay = -0; % This is the delay in msecs between the onset of
% the eye tracking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.

params.TrialInspectorFlag = 0; % If set, plot each time series

% Generally, only one or the other of the two harmonic flags should be set

params.HarmonicModelFlag = 1; % If set, model in the sparkline plots the harmonic
params.HarmonicTestFlag = 0; % If set, generate all response plots assuming harmonic stimulation

params.adapt_length = 60; % Duration in seconds of the initial background adaptation peiod

params.sampling_frequency = 20;  % The samples per second of our final data vectors
params.full_trial_length = 45;  % The length in seconds of the full trial
%   acquisition
params.final_trial_length = 40; % The length of each data trial after processing,
%   taken from the tail of each full data trial.
%    That is, we discard some data from the start of each trial.
params.minimum_length_trial = 42; % Any trial shorter than this will be discarded

params.StutterErrorFlag = 0;     % If these are data from the initial TTF4D
% experiment, use a modified version of the LeastSquaresSpectralFit routine
% that accounts for the tiny error in the stimulus sinusoidal modulation
params.spike_remover_params(1) = 10; % The window size to use for spike removal
params.spike_remover_params(2) = 0.5; % The max acceptable proportion change in the window
params.spike_remover_params(3) = 4; % The max acceptable change in SD units in the window
params.sgolay_span = 5;         % Parameters of the SGolay interpolation prior
params.sgolay_polynomial = 4;    %   to interpolation
params.BadPercentChangeThreshold = .40;       % Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is
% certainly noise
params.BadNanThreshold = 0.1; % If more than this proportion
% of a trial time series is composed of NaNs, discard the trial
params.SaveDataFlag = true;
params.SavePlotFlag = true;
params.RelabelDirections = true;
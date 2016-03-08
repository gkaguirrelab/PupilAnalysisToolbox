function params = PupilAnalysisToolbox_GetDefaultParams(varargin)
% params = PupilAnalysisToolbox_GetDefaultParams(varargin)
%
% Function that returns default parameters for the pupil analysis.
%
% 2/3/16    ms      Written.

% Parse the input
p = inputParser;

p.addOptional('SaveDataFlag', true, @islogical); % If set, save data
p.addOptional('SavePlotFlag', true, @islogical); % If set, save plots
p.addOptional('RelabelDirections', true, @islogical); % If set, relabel directions
p.addOptional('TrialInspectorFlag', false, @islogical); % If set, plot each time series
p.addOptional('HarmonicModelFlag', true, @islogical); % If set, model in the sparkline plots the harmonic
p.addOptional('HarmonicTestFlag', false, @islogical); % If set, generate all response plots assuming harmonic stimulation
p.addOptional('StutterErrorFlag', false, @islogical); % If these are data from the initial TTF4D
% experiment, use a modified version of the LeastSquaresSpectralFit routine
% that accounts for the tiny error in the stimulus sinusoidal modulation

% Duration definitions
p.addOptional('StimOnsetDelay', 0, @isscalar); % This is the delay in msecs between the onset of
% the eye tracking time stamps and the start of the stimulus modulation.
% This duration of data is chopped from the front of each trial prior to
% the analysis to allow the phase measurements to be accurate.
p.addOptional('adapt_length', 60, @isscalar); % Duration in seconds of the initial background adaptation period
p.addOptional('sampling_frequency', 20, @isscalar); % The samples per second of our final data vectors
p.addOptional('full_trial_length', 45, @isscalar); % The length in seconds of the full trial
p.addOptional('final_trial_length', 40, @isscalar); % The length of each data trial after processing,
% taken from the tail of each full data trial.
% That is, we discard some data from the start of each trial.
p.addOptional('minimum_length_trial', 35, @isscalar); % Any trial shorter than this will be discarded
p.addOptional('spike_remover_params', [8 0.5 5], @isscalar); % 3-element vector with the following elements:
% 1- The window size [in samples] to use for spike removal
% 2- The max acceptable proportion change in the window
% 3- The max acceptable change in SD units in the window

% Add parameters for the velocity-based blink detection algorithm
p.addOptional('VelocityBlinkDetectionFlag', true, @islogical);
p.addOptional('VelocityOnsetThreshold', -0.05, @isscalar);
p.addOptional('VelocityOffsetThreshold', 0.05, @isscalar);
p.addOptional('VelocitySearchWindowSize', 10, @isscalar);
p.addOptional('VelocityMarginWindowSize', 30, @isscalar);
p.addOptional('VelocitySmoothingParam', 5, @isscalar);

p.addOptional('sgolay_span', 20, @isscalar); % Parameters of the SGolay interpolation prior to interpolation
p.addOptional('sgolay_polynomial', 7, @isscalar);
p.addOptional('BadPercentChangeThreshold', 0.8, @isscalar); % Any time point with an absolute change from
% the mean of greater than this proportion will be NaN-ed, as it is certainly noise
p.addOptional('BadNaNThreshold', 0.1, @isscalar); % If more than this proportion
% of a trial time series is composed of NaNs, discard the trial
p.addOptional('yLim', 0.4, @isscalar); % Plotting range for y axis
p.addOptional('xLim', 35, @isscalar); % Plotting range for y axis

% Extract and assign the inputs
p.parse(varargin{:});
params = p.Results;
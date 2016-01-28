function [FitAmp,FitPhase,BestFit] = LeastSquaresSpectralFit_StutterError(data,frequency,sampling_frequency)
% Fits a sum of sine and cosine to normalized pupillary constriction data.
% Assumes that constriction corresponds to an increase in signal, therefore
% all data are multiplied by -1.
%
% The phase convention is that a sine wave has a phase of 0.
%
% 10/28/13  spitschan
% xx/xx/13  gka
%
% See http://www.mathworks.com/matlabcentral/answers/82576 for more details 
%
% On May 4, 2013, MS discovered that the TTF4D dataset had been collected
% with an imperfection in the creation of the sinusoidal modulation of the
% stimulus. Details are here:
%
%    https://cfn.upenn.edu/aguirre/wiki/private:melanopsin_pupillometry_protocols_data#temporal_transfer_function_in_four_directions_ttf4d
%
% This modified version of the LeastSquaresSpectralFit accounts for this
% error, and is to be used only for the analysis of the TTF4D data
% collected prior to May 1, 2013.
%
% - GKA
%

if (isempty(sampling_frequency))
    sampling_frequency = 20; % assume 20 samples per second if not specified
end

% Create one cycle of predicted modulation with 200 steps
  
num_steps_in_cycle=(1/frequency)*sampling_frequency;
num_cycles=ceil(length(data)/num_steps_in_cycle);

s_one_cycle=sin(linspace(0,2*pi,200));
s_resample_one_cycle=resample(s_one_cycle,num_steps_in_cycle,200);
s=repmat(s_resample_one_cycle,1,num_cycles);

c=fshift(s,num_steps_in_cycle/4);

% Trim from full cycles down to the actual data length

s=s(1:length(data));
c=c(1:length(data));


% To calculate phase properly, we multiple the vector by (-1). This
%  aligns the positive direction of change of the stimulus with a
%  positive change in the signal.

vec = data*(-1);

G=transpose([c;s]);
Betas=regress(vec',G);
FitAmp = sqrt(sum(Betas.^2));
FitPhase = atan2(Betas(1),Betas(2));

% The BestFit is multiplied by (-1) to return it to pupil contraction
% values

BestFit = ((c.*Betas(1))+(s.*Betas(2)))*(-1);

end

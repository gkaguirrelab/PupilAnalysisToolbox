function [FitAmp,FitPhase,BestFit] = LeastSquaresSpectralFit(data,frequency,sampling_frequency)
% Fits a sum of sine and cosine to normalized pupillary constriction data.
% Assumes that constriction corresponds to an increase in signal, therefore
% all data are multiplied by -1.
%
% The phase convention is that a sine wave has a phase of 0.
%
% See http://www.mathworks.com/matlabcentral/answers/82576 for more details 
%
% 10/28/13  spitschan
% xx/xx/13  gka

if (isempty(sampling_frequency))
    sampling_frequency = 20; % assume 20 samples per second if not specified
end

% To calculate phase properly, we multiple the vector by (-1). This
%  aligns the positive direction of change of the stimulus with a
%  positive change in the signal.
vec = data*(-1);
c=cos(linspace(0,2*pi*length(vec)/sampling_frequency*frequency,length(vec)));
s=sin(linspace(0,2*pi*length(vec)/sampling_frequency*frequency,length(vec)));
G=transpose([c;s]);
Betas=regress(vec',G);
FitAmp = sqrt(sum(Betas.^2));
FitPhase = atan2(Betas(1),Betas(2));

% The BestFit is multiplied by (-1) to return it to pupil contraction
% values

BestFit = ((c.*Betas(1))+(s.*Betas(2)))*(-1);

end

function iy = SGolaySmooth(times,values,fitspan,poly,sampling_frequency,full_trial_length)

if (isempty(sampling_frequency))
    sampling_frequency = 20; % assume 20 samples per second if not specified
end

x = times/1000;           % convert measurement times to secs
y = values;               % place values in a new variable


% To avoid a spike in fitting at the vector onset, add a set of initial
% values which are equal to the first recorded value.

if (x(1)>0.1)
    SamplesToPad=round(x(1)-0.05)/0.05;
    InitialPadX=linspace(0,x(1)-0.05,SamplesToPad);
    x=[InitialPadX,x];
    y=[ones(1,SamplesToPad).*y(1),y];
end

sy = smooth(x,y,fitspan,'sgolay',poly);  % smoothed data
ix = linspace(0,full_trial_length-(1/sampling_frequency),full_trial_length*sampling_frequency);
iy = interp1(x,sy,ix,'spline');  % interpolated data to even sample times, starting from zero
end
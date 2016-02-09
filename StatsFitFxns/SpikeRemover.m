function [cleaned, indx] = SpikeRemover(data,spike_remover_params)
% [cleaned, indx] = SpikeRemover(data,spike_remover_params)
%
% This routine takes in three arguments in remover_params:
%   1- The window size to use for spike removal
%   2- The max acceptable proportion change in the window
%   3- The max acceptable change in SD units in the window
% This routine is pretty bone-headed. It sets to NaN any stretch of WINDOW
% measures from the data within which the data contains a max and min that
% are more than a X% signal change or Y SD of the entire time series.
% The WINDOW before and after are also set to NaN

window = spike_remover_params(1); % the size of the window to sweep
spike_thresh_percen = spike_remover_params(2); % Remove points of greater than X% change
spike_thresh_sd = spike_remover_params(3); % Remove points of greater than X SD change

start_spot=floor(window*1.5)+1;
end_spot=length(data)-floor(window*1.5)-1;

cleaned = data;
std_data = std(data);
indx = [];
for q = start_spot:end_spot
    if ( abs(nanmax(data(q-floor(window/2):q+floor(window/2))) - nanmin(data(q-floor(window/2):q+floor(window/2)))) > spike_thresh_percen)
        cleaned(q-floor(window*1.5):q+floor(window*1.5))=NaN;
        indx = [indx q-floor(window*1.5):q+floor(window*1.5)];
    end % if statement for perecentage change
    if ( abs(nanmax(data(q-floor(window/2):q+floor(window/2))) - nanmin(data(q-floor(window/2):q+floor(window/2)))) > spike_thresh_sd*std_data)
        cleaned(q-floor(window*1.5):q+floor(window*1.5))=NaN;
        indx = [indx q-floor(window*1.5):q+floor(window*1.5)];
    end % if statement for SD change
    
   % 1) what's and where the min value in time series
   % 2) what's the max value before the point, and after the point
   % 3) calculate diff between each max and min
   % 4) take smaller of the two diffs
   % 5) is that difference greater than the threshold?
   % 6) code: flag ? bidirectional vs. unidirectional
    
end % move window through data

end % function
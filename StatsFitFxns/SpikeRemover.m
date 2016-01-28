function [cleaned, indx] = SpikeRemover(data,spike_remover_params)

% This routine is pretty bone-headed. It sets to Nan any stretch of WINDOW
% measures from the data within which the data contains a max and min that
% are more than a 20% signal change or 3 SD of the entire time series.
% The WINDOW before and after are also set to NAN

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
    
    
    
end % move window through data

end % function
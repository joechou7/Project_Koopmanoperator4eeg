function [winstart,nw,mod_end,interval] = slicewin(n_ts,winsize,overlap)
% This function pre-processes time series by slicing signal into windows
% and create indicators for change points
%   Input:
% n_ts - the length of time series
% winsize - the window size
% overlap - the overlapping rate in [0,1]
%   Output:
% winstart - a vector to indicate the start point of each window
% nw - number of windows
% mod_end - the modified end point of time series (the last points which are not sufficient to form a window
% will be abandoned by default)

overlap_zone = floor(winsize*overlap);
nw = floor((n_ts - winsize)/(winsize - overlap_zone))+1; % The number of windows
interval = winsize - overlap_zone;
winstart = 1:interval:(n_ts - winsize + 1);
mod_end = winstart(end) + winsize - 1;
end
function res = label_seizure(start,ends,windowsize,thre_precent)
%LABEL_SEIZURE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    thre_precent = 0.25;
end


whole_s = start/windowsize;
rest_s = whole_s-floor(whole_s);
if rest_s > 1-thre_precent
    start = floor(whole_s) + 1;
else
    start = floor(whole_s);
end

whole_e = ends/windowsize;
rest_e = whole_e - floor(whole_e);
if rest_e < thre_precent
    ends = floor(whole_e);
else
    ends = floor(whole_e)+1;
end
res = start:ends;

end


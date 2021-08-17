function D = system_recon(attractor, target, ev, method)
%SYSTEM_RECON Summary of this function goes here
%   Detailed explanation goes here
% target: t * n
% D: t * 1
if nargin < 4
    method = 'wteud';
end
if strcmp(method, 'wteud')
    D = sqrt(((target-attractor).^2) * ev);
elseif strcmp(method, 'max')
    D = max(abs(target-attractor),[],2);
elseif strcmp(method, 'sin')
    
end

end


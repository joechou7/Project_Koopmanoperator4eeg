function [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(data,dt, stackmax, rmax)
%HANKLE_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here

% H = zeros(stackmax,size(data,1)-stackmax);
% parfor k=1:stackmax
%     H(k,:) = data(k:end-stackmax-1+k,1);
% end

index1 = 1:stackmax;
index2 = stackmax:size(data,1)-1;
c = data(index1,1).'; r = data(index2,1);
H = hankel(c,r);

[U,S,V] = svd(H,'econ');
sigs = diag(S);
beta = size(H,1)/size(H,2);
tspan = dt:dt:dt*length(data);
x = [];
y = [];
t = [];
r = rmax;

end


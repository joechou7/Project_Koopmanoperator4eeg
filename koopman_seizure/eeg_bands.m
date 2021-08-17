function [cleandata,filtdelta,filttheta,filtalpha,filtbeta,filtgamma] = eeg_bands(rawdata,order)
%EEG_BANDS Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    order = 3;
end
Fs = 256;
filtdata = eegfilt(rawdata,Fs,0.5,0,length(rawdata),order);
cleandata = eegfilt(filtdata,Fs,0,50,length(rawdata),order);
filtdelta = eegfilt(cleandata,Fs,0.5,4,length(rawdata),order);
filttheta = eegfilt(cleandata,Fs,4,8,length(rawdata),order);
filtalpha = eegfilt(cleandata,Fs,8,13,length(rawdata),order);
filtbeta = eegfilt(cleandata,Fs,13,30,length(rawdata),order);
filtgamma = eegfilt(cleandata,Fs,30,48,length(rawdata),order);
end


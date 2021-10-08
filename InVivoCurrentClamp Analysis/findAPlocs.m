function [pkLocs, pkLocsByTrace,XAPs, Xsub, APThreshold] = findAPlocs(X, samplingFreq, minPeakHeight, peakDetectionHeightDiscountFactor, mpw)
%findAPlocs - finds AP threshold and locations, and parses APs from the
%underlying subthreshold potential
%   X can be a nFrames x nFiles matrix;  retruns pkLocs assuming
%   concatenation of all files to a single row vector (else return
%   pkLocsByTrace), APs holds a matrix of extracted APs, Xsub returns
%   concatenated input traces with APs removed, XAP is a binary vector of
%   length X with AP locs
%% paramaters
mpwFrame = round(mpw/1000*samplingFreq); %convert min peak width for AP detection from ms to frames
[~, nFiles] = size(X);
totalElements = numel(X);


% diagnosticPlotFindPeaks=1;
%% find APs
Xnorm = X-quantile(X,0.05); %normalize traces
concatX = reshape(Xnorm,1,[]); %concat traces
maxPeak = max(concatX); %find max peak and find APs of similar size
pkLocsByTrace = [];
[~,pkLocs] = findpeaks(concatX, 'MinPeakHeight', max([minPeakHeight maxPeak*peakDetectionHeightDiscountFactor]), 'MinPeakWidth', mpwFrame);
nPks = length(pkLocs);
%% generate AP vector
XAPs = zeros(1,totalElements); %holds APs as 0 and 1s
XAPs(pkLocs)=1;
Xsub = concatX; %holds subthreshold mV with AP removed
APThreshold = nan(nPks,1); %will hold threshold value for AP
thresholdLocs = nan(nPks,1); %hold location of where threshold is tripped
if nPks>1
[APThreshold, thresholdLocs] = getAPThreshold(Xsub,pkLocs,samplingFreq);
end

%% remove spikes from subthreshold recording (adapted from Azouz & Gray, 1999)
%removes a region around spike and replaces with NaN
for iPk = 1:nPks
    try %will nto work work if thresholdLoc is NaN
    startPt = thresholdLocs(iPk); %start from threshold
    endPt =find(Xsub(pkLocs(iPk):min([totalElements pkLocs(iPk)+maxLookAheadFrames]))<Xsub(thresholdLocs(iPk)),1,'first')-1; %where ap returns below threshold
        if isempty(endPt)
        endPt = pkLocs(iPk)+maxLookAheadFrames;
        else
        endPt = pkLocs(iPk) + endPt;
        end
    Xsub(startPt:endPt) = NaN;
    end
end

Xsub = reshape(Xsub,[],nFiles);
XAPs = reshape(XAPs,[],nFiles);

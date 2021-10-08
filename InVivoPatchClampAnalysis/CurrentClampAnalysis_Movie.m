%% CurrentClampAnalysis_Movie
%analyzes current clamp response to natural movie, mainly AP spiking (meanAPFreq)
%% Paramaters
close all;
clear all;

durationOfDrift_sec =2.5; %duration of drifting grating
durationOfGrey_sec =0.5;%1;
nAngles = 12; %number of angles including drift direction
angles = linspace(0,330,nAngles);
DesiredSamplingFreq = 20000; %or desired in the case of downsampling

%paramaters for AP detection 
mpw = 0.25; %for peak detection in msec
peakDetectionHeightDiscountFactorDefault = 0.75; %how much off of max height of tallest spike can other spikes be detected
minPeakDetectionHeightDefault = 15;
%% Current Clamp Analysis
%get traces from pre-defined RecordingFilePaths
%for each Recording do:

%initialize variables
RMP = nan(nRecordings,1); %contains mean resting membrane potential from natural movie exps
VSpikeThreshold = nan(nRecordings,1); %contains mean spike threshold from  natural movie exps
meanAPFreq = nan(nRecordings,1); %contains mean APfreq
meanSub = nan(nRecordings,1); %contains mean depol

for iRecording = 1:nRecordings
[X,samplingFreq,abfTimeStamp] = getABFfiles(RecordingFilePaths,DesiredSamplingFreq);
%X holds primary channel (current or voltage) 
[meanSub(iRecording), meanAPFreq(iRecording) ,VspikeThresholdMean(iRecording)] = visMovieEphysAnalysis(X,samplingFreq,peakDetectionHeightDiscountFactor,minPeakHeight,mpw)
RMP(iRecording) = nanmedian(quantile(X,0.05,1));
end

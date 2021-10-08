%% VoltageClampAnalysis_Grating
%analyzes grating responses in voltage clamp, mainly the mean, std, cv of
%responses, and orientation tuning width

%% Paramaters
close all;
clear all;

durationOfDrift_sec =2.5; %duration of drifting grating
durationOfGrey_sec =0.5;%1;
nAngles = 12; %number of angles including drift direction
angles = linspace(0,330,nAngles);
DesiredSamplingFreq = 20000; %or desired in the case of downsampling

%initialize data for storage
meanResponseTuningCurveCentered = nan(nRecordings,nAngles/2); %contains subthreshold tuning curve
stdResponseTuningCurveCentered = nan(nRecordings,nAngles/2); %contains trial-to-trial stdev of subthreshold tuning curve
CVResponseTuningCurveCentered = nan(nRecordings,nAngles/2); %contains trial-to-trial CV of subthreshold tuning curve
tuningWidth = nan(nRecordings,1);  %contains sigma for gaussian fit of subthrehsold tuning curve, sigma = 1 stdev

%% Extract Traces
%get traces from pre-defined RecordingFilePaths
%X holds primary channel (current or voltage) 
%for each Recording do:
[X,samplingFreq,abfTimeStamp] = Zahid_getABFfiles(RecordingFilePaths,DesiredSamplingFreq);
%% Get mean, std,  CV of current responses per angle for each cell, and tuningWidth

for iRecording = 1:nRecordings
           [meanResponseTuningCurveCentered(iRecording,:),stdResponseTuningCurveCentered(iRecording,:),CVResponseTuningCurveCentered(iRecording,:),tuningWidth(iRecording)] = visStimEphysAnalysis_vclamp(X,samplingFreq,angles,angleTimingsSorted,durationOfDrift_sec, durationOfGrey_sec);
end


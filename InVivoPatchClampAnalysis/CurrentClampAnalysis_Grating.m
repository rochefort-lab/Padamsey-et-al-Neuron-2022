%% CurrentClampAnalysis_Grating
%analyzes grating responses in current clamp, mainly the mean, std, cv of
%responses, and orientation tuningWidth
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
    
%set up matrix for data storage
Analysis = []; %structure to store analysis
Analysis.Stim.meanTuningCurve.subThresh = nan(nRecordings,nAngles/2); %contains subthreshold tuning curve
Analysis.Stim.meanTuningCurve.APFreq = nan(nRecordings,nAngles/2); %contains AP tuning curve
Analysis.Stim.stdTuningCurve.subThresh = nan(nRecordings,nAngles/2); %contains trial-to-trial stdev of subthreshold tuning curve
Analysis.Stim.stdTuningCurve.APFreq = nan(nRecordings,nAngles/2); %contains trial-to-trial stdev of AP tuning curve
Analysis.Stim.CVTuningCurve.subThresh = nan(nRecordings,nAngles/2); %contains trial-to-trial CV of subthreshold tuning curve
Analysis.Stim.CVTuningCurve.APFreq = nan(nRecordings,nAngles/2);  %contains trial-to-trial CV of AP tuning curve
Analysis.Stim.meanDistToSpikeTuningCurve = nan(nRecordings,nAngles/2); %contains meanDistanceToSpike (distance from peak response to spike threshold)
Analysis.Stim.stdDistToSpikeTuningCurve = nan(nRecordings,nAngles/2); %contains trial-to-trial stdev of meanDistanceToSpike (distance from peak response to spike threshold)
Analysis.Stim.CVDistToSpikeTuningCurve = nan(nRecordings,nAngles/2); %contains trial-to-trial CV of meanDistanceToSpike (distance from peak response to spike threshold)
Analysis.Stim.nStimFiles = nan(nRecordings,1); %number of stim files
Analysis.Stim.OSIGaussian.APFreq = nan(nRecordings,1); %contains sigma for gaussian fit of AP tuning curve, sigma = 1 stdev
Analysis.Stim.OSIGaussian.subThresh = nan(nRecordings,1);  %contains sigma for gaussian fit of subthrehsold tuning curve, sigma = 1 stdev
Analysis.Stim.VSpikeThreshold = nan(nRecordings,1); %contains mean spike threhsold from orientation grating exps
Analysis.RMP = nan(nRecordings,1); %contains mean resting membrane potential from orientation grating exps
Analysis.Movie.RMP = nan(nRecordings,1); %contains mean resting membrane potential from natural movie exps
Analysis.Movie.VSpikeThreshold = nan(nRecordings,1); %contains mean spike threshold from  natural movie exps
Analysis.Stim.trialSubDistAPData = cell(nRecordings,1); %holds subthrehold depol (mV) and spike (binary 0,1) values

SubOrAPlabels = ["subThresh" "APFreq"];
nSubOrAPlables = numel(SubOrAPlabels);
%% Extract Traces
%get traces from pre-defined RecordingFilePaths
%X holds primary channel (current or voltage) 
%for each Recording do:
[X,samplingFreq,abfTimeStamp] = Zahid_getABFfiles(RecordingFilePaths,DesiredSamplingFreq);
%% Get mean, std, and CV of subthreshold and AP responses per angle for each cell

%for each recording do the following:

for iRecording = 1:nRecordings
            [meanResponseTuningCurve_recording, stdResponseTuningCurve_recording, CVResponseTuningCurve_recording,...
                meanDistToSpikeTuningCurve_recording, stdDistToSpikeTuningCurve_recording, CVDistToSpikeTuningCurve_recording,...
                Analysis.Stim.VSpikeThreshold(iRecording), Analysis.Stim.trialSubDistAPData{iRecording},Analysis.Stim.GaussianR.APFreq{iRecording},...
                Analysis.Stim.GaussianR.subThresh{iRecording}]
                = visStimEphysAnalysis(X,samplingFreq, angles, angleTimingsSorted, durationOfDrift_sec; durationOfGrey_sec,peakDetectionHeightDiscountFactor(iRecording),minPeakDetectionHeight(iRecording),mpw);
                   
            Analysis.RMP(iRecording) = nanmedian(quantile(X(1:30000,1),0.05));
            Analysis.Stim.meanTuningCurve.subThresh(iRecording,:) = meanResponseTuningCurve_recording(1,:);
            Analysis.Stim.meanTuningCurve.APFreq(iRecording,:) = meanResponseTuningCurve_recording(2,:);
            Analysis.Stim.stdTuningCurve.subThresh(iRecording,:) = stdResponseTuningCurve_recording(1,:);
            Analysis.Stim.stdTuningCurve.APFreq(iRecording,:) = stdResponseTuningCurve_recording(2,:);
            Analysis.Stim.CVTuningCurve.subThresh(iRecording,:) = CVResponseTuningCurve_recording(1,:);
            Analysis.Stim.CVTuningCurve.APFreq(iRecording,:) = CVResponseTuningCurve_recording(2,:);
            Analysis.Stim.meanDistToSpikeTuningCurve(iRecording,:) = meanDistToSpikeTuningCurve_recording;
            Analysis.Stim.stdDistToSpikeTuningCurve(iRecording,:) = stdDistToSpikeTuningCurve_recording;
            Analysis.Stim.CVDistToSpikeTuningCurve(iRecording,:) = CVDistToSpikeTuningCurve_recording;
            Analysis.Stim.nStimFiles(iRecording,1) = nAbfStimFiles;
            Analysis.Stim.meanTuningCurveNorm.APFreq=Analysis.Stim.meanTuningCurve.APFreq./max(Analysis.Stim.meanTuningCurve.APFreq,[],2);
            Analysis.Stim.meanTuningCurveNorm.subThresh=Analysis.Stim.meanTuningCurve.subThresh./max(Analysis.Stim.meanTuningCurve.subThresh,[],2);
end
%% Define spikers
%a spiker is a cell that spikes at least once in response to visual
%stimulation, only spikers are used for AP related analyses
spikers = max(Analysis.Stim.meanTuningCurve.APFreq,[],2)>0; %only imposed for spike output tuning width, 
%% Orientaiton tuning widths (bassed on Gaussian fits)
%tuning width = 1 stdev of tuning curve based on Gaussian fit
for iRecording = 1:nRecordings
    for iSubOrAP = 1:nSubOrAPlables %fit gaussian to subthreshold or ap curves
        X= Analysis.Stim.GaussianR.(SubOrAPlabels{iSubOrAP}){iRecording};        
            if iSubOrAP==2 %ie AP tuning fits
                RTSC = GaussianFit_Orientation(1,X,0); %RTSC is 4 paramaters of gaussian, 3rd paramater is tuning width
            else %ie for subthreshold tuning fits
                RTSC = GaussianFit_Orientation(1,X,0,0, [30 45 90 135 180 225 270 315 360]); %last variable holds set of sigma values for initialization
            end
            Analysis.Stim.OSIGaussian.(SubOrAPlabels{iSubOrAP})(iRecording) = RTSC(3);
    end
end
Analysis.Stim.OSIGaussian.(SubOrAPlabels{2})(~spikers) = NaN; %omit AP tuning values for cells that did not spike
%% Calculate I/O curve (distToSpike vs probability/frequency of spike)
%this examines full ephys trace and partions recording to bins to determine
%average membrane depolarization and average spike probability/frequency
%bins defined by:
integrationTimeForAP = 50; %in ms
integrationFramesForAP = convertMillisecondsToFrames(integrationTimeForAP,samplingFreq);
binarize = 1; %if 1, I/O is depol vs p(spike) else it is depol vs spike freq

%once bins of subthreshold mv and AP are made, we need to place subthrehsold
%mv in to bins in order to calculate the probability of observing a spike
%within a range of depolarization 

%paramaters for this binning: 
nHistBins = 10;
edges = linspace(0,1,nHistBins+1);
pSpike = nan(numel(spikers),nHistBins); % initialize - will hold I/O curve (pSpike vs subthresold depol for each file (ie spiker, which is a cell that spikes in response to vis stim))


for iRecordings = find(spikers)'
    APThreshold = nanmean(Analysis.Stim.APThreshold{iRecording}); %AP threshold for the recording
    X=Analysis.Stim.trialSubDistAPData{iRecording}; %trialSubDistAPData = [Xsub,XAP]
    
    %linearize
    nFiles = size(X,2)/2;
    X(:,1:nFiles)=X(:,1:nFiles)./APThreshold; %standardize distance to spike, makes comparison between animals and fd vs con group fair  
    Xsub = reshape(X(:,1:nFiles),[],1);%./nanmean(Analysis.Stim.meanAPThresholdTuningCurve(iRecording,:));
    XAP = reshape(X(:,nFiles+1:end),[],1);

    meanDepol = movmin(Xsub,integrationFramesForAP); %mean subthreshold depol within time window
    XAP = movmean(XAP,integrationFramesForAP); %check to see if depol generated spike
      
    if binarize
        XAP(XAP>0) = 1; %binarize within a bin we just care if there was at least one spike for calculating p(spike) else for spike freq binarize should be 0
    end
    
    
    [NperBin,~,binID] = histcounts(meanDepol,edges); %use histcounts to bin subdepol bins based on value defined by edges, and calculate p(spike)
    for ibin = 1:nHistBins
        val = XAP(find(binID==ibin));
        pSpike(iRecording,ibin) = sum(val)/(numel(val)); %p(spike)
    end
end
%% Calculation of trial-to-trial variability (CV)
%takes CV of subthrehsold variability (subthresholdCV) where the mean response is more than
%10mV below spike threshold, this is to avoid a ceiling effect (ie
%variability is reduced for strong responses as they always spike, meaning
%their max subthreshold depol cannot exceed AP threshold)

subthresholdCV = Analysis.Stim.CVTuningCurve.subThresh.*(Analysis.Stim.meanDistToSpikeTuningCurve<-10); %only take responses >10mV from threshold
subthresholdCV(subthresholdCV==0) = NaN;
subthresholdCV = nanmean(subthresholdCV,2); %for each cell, average CV across orientations

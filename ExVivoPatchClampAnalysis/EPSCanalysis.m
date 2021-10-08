%% EPSCanalysis
%calculates EPSC amplitude, PPR, invCV2 on ex vivo EPSC data
%% Paramaters
close all;
clear all;

samplingFreq = 20000;
stimArtefactLocsThreePulse_ms = ([1158.33;5163.03;5213.03]); %3 stim pulses; 1st to assess i/o, 2nd and 3rd pulses assess PPR
PPR = nan(nRecordings,1); %paired pulse ratio
invCV2 = nan(nRecordings,1); %inverse CV^2 measure
meanResponseAmplitude = nan(nRecordings,1); %mean response for making i/o curve 

%% get EPSC amplitudes
for iRecording = 1:nRecordings
            X = squeeze(abf2load(char(RecordingFilePath(iRecording))));
            %send only primary Channel
            X = squeeze(X(:,1,:)); %returns data for primary channel of recording. Each recording as several sweeps, each of which has 3 stim pulses; the last two pulses are delivered as paired pulses for PPR analysis
            [responseAmplitude] = getResponseAmplitudes(X, samplingFreq, stimArtefactLocsThreePulse_ms); %responseAmplitude for trace x pulse
            PPR(iRecording) = nanmean(responseAmplitude(:,3),1)/nanmean(responseAmplitude(:,2),1);
            invCV2(iRecording) = (nanmean(responseAmplitude(:,1),1)/nanstd(responseAmplitude(:,1),[],1))^2;
            meanResponseAmplitude(iRecording) = nanmean(responseAmplitude(:,1),1);
end     

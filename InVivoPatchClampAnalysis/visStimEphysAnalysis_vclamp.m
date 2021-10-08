function [meanResponseTuningCurveCentered,stdResponseTuningCurveCentered,CVResponseTuningCurveCentered,tuningWidth] = visStimEphysAnalysis_vclamp(X,samplingFreq,angles,angleTimingsSorted,durationOfDrift_sec, durationOfGrey_sec)
%% paramaters
nAngles = length(angles); %includes directions

%rectify and subtract baseline
X=sign(sum(X,1)).*X; %allows same function to be used for inhibition/excitation
X=X-quantile(X,0.05);

[~, nFiles] = size(X);
baselineDurationFrames = round(durationOfGrey_sec*samplingFreq);
driftDurationFrames = round(durationOfDrift_sec*samplingFreq);
angleFramingsSorted = round(angleTimingsSorted*samplingFreq);

%% Segment Files based on angleTimings
Xsorted = nan(nFiles,nAngles,driftDurationFrames+baselineDurationFrames); %files x angles x time (stim response period)
for iFile = 1:nFiles
    [~,idx] = min(angleFramingsSorted(:,iFile));
    angleFramingsSorted(idx,iFile)=NaN;
    for iAngle = 1:nAngles
     Xsorted(iFile,iAngle,:) = X(angleFramingsSorted(iAngle,iFile)-baselineDurationFrames:angleFramingsSorted(iAngle,iFile)+driftDurationFrames-1,iFile);  
    end
end
%% calculate mean and stdev tuning curves and OSI

quantileResponse = 0.5; %take median
responses = quantile(Xsorted(:,:,baselineDurationFrames:end),quantileResponse,3)-quantile(Xsorted,quantileResponse,3); %take median response within response window
meanResponseTuningCurve = nanmedian(responses,1); %mean across trials
meanResponseTuningCurve = (meanResponseTuningCurve(1:6)+meanResponseTuningCurve(7:12))/2; %fold responses (ie average directions of the same orientaiton)

%calc stdev
stdResponseTuningCurve = nanstd(responses,[],1);
stdResponseTuningCurve = (stdResponseTuningCurve(1:6)+stdResponseTuningCurve(7:12))/2; %fold responses

%center tuning curve
[~, centeringIdxs] = CenterTuningCurve(meanResponseTuningCurve);

meanResponseTuningCurveCentered = nan(nAngles/2);
stdResponseTuningCurveCentered = nan(nAngles/2);
for iAngle=1:nAngles/2
meanResponseTuningCurveCentered(centeringIdxs(iAngle)) = meanResponseTuningCurve(iAngle);
stdResponseTuningCurveCentered(centeringIdxs(iAngle)) = stdResponseTuningCurve(iAngle);
end

CVResponseTuningCurveCentered = (stdResponseTuningCurveCentered./meanResponseTuningCurveCentered);
CVResponseTuningCurveCentered(meanResponseTuningCurveCentered<10)=NaN; %if mean <10 pA ignore else noisy

%determine tuning width 
RTSC = GaussianFit_Orientation(1,responses,0,0, [30 45 90 135 150 180]); %last variable holds set of sigma values for initialization
tuningWidth = RTSC(3);%tuning width is 3rd variable of output
end

%%

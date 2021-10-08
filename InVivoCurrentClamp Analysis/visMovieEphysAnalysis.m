function [meanSub meanAPFreq ,VspikeThresholdMean] = visMovieEphysAnalysis(X,samplingFreq,peakDetectionHeightDiscountFactor,minPeakHeight,mpw)
%returns means subthreshold depolarization (meanSub), firing frequency
%(meanAPFreq) and spike threshold VspikeThresholdMean

%% paramaters
    [nTotalFrames nFiles] = size(X);
   
     %get APs
    [pkLocs, ~,  XAPs, Xsub, APThreshold] = findAPlocs(X, samplingFreq, minPeakHeight, peakDetectionHeightDiscountFactor, mpw);
    nPks = length(pkLocs); %number of APs in file
    VspikeThresholdMean = nanmedian(APThreshold); %mean spike threshold

    % calcmean subthreshold depol
    %rectify and subtract baseline
    X=sign(sum(X,1)).*X; %allows same function to be used for inhibition/excitation
    X=X-quantile(X,0.05); %remove baseline
    reshapeX = reshape(Xnorm,1,[]);
     meanSub = nanmean(nanmean(Xsub));%meand depolarization for duration

    % calc spike freq
    recordingDuration =sum(~isnan(reshapeX))/samplingFreq;
    meanAPFreq = nPks/recordingDuration;
    
    
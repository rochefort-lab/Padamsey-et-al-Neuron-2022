function [APThreshold, APThresholdLoc] = getAPThreshold(Xsub,pkLocs,samplingFreq)
    %finds AP threhsold based on kink method (i.e. maximizing 2nd derivative)    
    lookBackTime_ms = 5; %lookback time from peak)
    lookBackFrames = round(lookBackTime_ms*samplingFreq/1000);
    
    nPks = numel(pkLocs);
    APThresholdLoc = nan(nPks,1);
    for iPk=1:nPks
    startPt = max([1 pkLocs(iPk)-lookBackFrames]);
    
    %finds location within lookbacktime which maximizes 2nd derivative,
    %this location is the threshold
    preAP = Xsub(startPt:pkLocs(iPk)); 
    diffPreAP = preAP(2:end)-preAP(1:end-1);
    diffdiffPreAP = diffPreAP(2:end)-diffPreAP(1:end-1);
    [~,maxChangeLoc] = max(diffdiffPreAP);
    maxChangeLoc = max([1 maxChangeLoc-2]); %correction shift;  
    APThresholdLoc(iPk) = startPt+maxChangeLoc-1;
    end
    
    APThreshold = Xsub(APThresholdLoc);

end


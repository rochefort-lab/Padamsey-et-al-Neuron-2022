function [responseAmplitude] = getResponseAmplitudes(X,samplingFreq, artefactLocs_ms)
%EVOKEDSTIMANALYSIS S
%get response amplitude of EPSCs
%% Paramaters

AmplitudeWindow = 15;%15; %ms window for defining response
AmplitudeWindowFrames = round(AmplitudeWindow/1000*samplingFreq);

lookBackWindow = 3; %ms duraiton to lookback from stim artefact for Baseline
lookBackWindowFrames = round(lookBackWindow/1000*samplingFreq);

BaselineWindow = 5; %ms window for calculating baseline
BaselineWindowFrames= round(BaselineWindow/1000*samplingFreq);

artefactLocsFrames = round((artefactLocs_ms)/1000*samplingFreq);
stimLocsFrames   = artefactLocsFrames; %use user input of artefact locations as stim location

nTraces = size(X,2);


%% EPSC amplitude
response = nan(nTraces,3); %traces,pulses
baseline = nan(nTraces,3);
targetLocsFrames = [stimLocsFrames,nullLocsFrames'];
    for iPulse = 1:3 %3 stim pulses are present per trial pulse 2 and 3 are paired pulses for assessing PPR.  
        for iTrace = 1:nTraces %go oneby one for response amplitude in case ap
            a = X(targetLocsFrames(iPulse,1):targetLocsFrames(iPulse,1)+AmplitudeWindowFrames,iTrace);
            [response(iTrace,iPulse)] =max(a); %chose peak response within window
        end
           baseline(:,iPulse) = max(X(targetLocsFrames(iPulse,iSignal)-lookBackWindowFrames-BaselineWindowFrames:targetLocsFrames(iPulse,iSignal)-lookBackWindowFrames,:)); %choose peak baseline response
    end

responseAmplitude = response-baseline; %traces x pulses

end


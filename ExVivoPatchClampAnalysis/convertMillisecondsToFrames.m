function [frames] = convertMillisecondsToFrames(time,samplingFreq )
%CONVERTTIMETOFRAMES 
frames = round(time*samplingFreq/1000);
end


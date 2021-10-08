%% VoltageClampAnalysis_Movie
%analyzes movie responses in current clamp, namley meanCurrent and
%meanNaFlux which are used for ATP calculations


%get traces from pre-defined RecordingFilePaths
%for each Recording do:

%initialize variables
meanCurrent = nan(nRecordings,1); %contains mean depol

for iRecording = 1:nRecordings
[X,samplingFreq,abfTimeStamp] = getABFfiles(RecordingFilePaths,DesiredSamplingFreq);
%X holds primary channel (current or voltage) 
meanCurrent(iRecording) = nanmean(nanmean(X));

% for excitatory current:
meanNaFlux = meanCurrent*1.42; %as per Harris et al., 2015 (Current Biology) x 1.42 accounts for Na/K overlap through AMPAR
end

%% Mini Analysis
% Analyzes ex vivo minis, getting amplitude (miniAmplitdue), frequency
% (miniFreq) and waveform (traceStore)
%% Parameters
close all;
clear all;
mpp = 1;%4; %min peak prominence for matlab peak detection
mpw = 1; %min peak width for matlab peak detection
mpd = 1; %min peak distance
min_mini_height = 1;%5; %min height filter applied at end, set to 1 if overrie
max_mini_height = 50; %max height filter applied at end
lookBackSize = 16; %how far to look back from peak to calculate baseline (will be multiplied for 1.5x for buffer)
lookFwdSize = 150; %how far to look for end (will be multiplied by 1.5x for buffer)
thresholdTemplateCorr = 0.6;
totalRecordingLengthinSeconds = 50*6;
sigma = 3; %no of thersholds above noise
samplingFreq = 20000;
%% Find and extract minis

%X = ephys trace loaded from metadata using abf2load
%first filter to remove hairness which inteferes with peak detection
lowpass_filter = designfilt('lowpassfir','FilterOrder',10,'CutoffFrequency',200, ...
            'DesignMethod','window','Window',{@kaiser,3},'SampleRate',20000);
        X_signal_filt = filter(lowpass_filter,X);

%find peaks
[~, signalLocs] = findpeaks(X_signal_filt, 'MinPeakProminence', mpp, 'MinPeakWidth', mpw, 'MinPeakDistance', mpd);
nLocs = numel(signalLocs);
miniFreq = nLocs/totalRecordingLengthinSeconds;

%Extraction of mini and calculation of peak amplitude (works better than matlab)
miniAmplitude = nan(nLocs); %mini Amplitude
baselineError = nan(nLocs);
        
traceLength = lookBackSize*1.5 + lookFwdSize*1.5 + 1;
traceStore = nan(traceLength,nLocs,2);

                X_signal = X_signal_filt;
                X_raw = X;
                peakLocs = signalLocs;
             
            parfor iLoc = 1:nLocs
                
                             startPt = peakLocs(iLoc)-1.5*lookBackSize;
                    endPt = peakLocs(iLoc)+1.5*lookFwdSize;
                    trace = X_signal(startPt:endPt);
                    baseline = (min(X_signal(startPt:startPt+lookBackSize-1,1)));
                    detectedPeakProminence = X_signal(peakLocs(iLoc))-baseline;
                    
                    
                    %check for best peak
                    [potentialEarlyPeak,potentialEarlyPeakLoc] =(max(X_signal(startPt:startPt+lookBackSize-1,1)));
                    potentialEarlyBaseline = min(X_signal(startPt+potentialEarlyPeakLoc-1-1.5*lookBackSize:startPt+potentialEarlyPeakLoc-1+lookBackSize-1));
                    potentialEarlyPeakProminence = potentialEarlyPeak-potentialEarlyBaseline;
                    
                                        
                    %update locaiton
                    if potentialEarlyPeakProminence>detectedPeakProminence && abs(peakLocs(iLoc) - (startPt+potentialEarlyPeakLoc-1))<5
                        miniAmplitude(iLoc) = potentialEarlyPeakProminence;
                        peakLocs(iLoc) = startPt+potentialEarlyPeakLoc-1;
                        startPt = peakLocs(iLoc)-1.5*lookBackSize;
                        endPt = peakLocs(iLoc)+1.5*lookFwdSize;
                        trace = X_signal(startPt:endPt);
                        potentialEarlyBaseline = baseline;
                        signalLocs(iLoc) =  peakLocs(iLoc);
                    else
                        miniAmplitude(iLoc) = detectedPeakProminence; %peak amplitudes 
                        
                    end

                   traceStore(:,iLoc) = trace; %store trace for MeanVarianceAnalysis, conducted as per described in Padamsey et al., 2021 (Neuron)

            end
     
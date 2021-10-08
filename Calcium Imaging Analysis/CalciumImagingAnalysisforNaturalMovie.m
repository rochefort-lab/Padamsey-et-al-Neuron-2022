%% CalciumImagingAnalysisforNaturalMovie
%analyzes calcium imaging data for natural movies, mainly used for decoding
%information from inferred spikes (stored in decoderAccuarcy)

%% Paramaters
close all;
clear all;
imagingFreq = 20;

%frames sizes
binSize = 20; %default 20 frames, which is 1s; how to cut movie in frames *make multiples of movieSize factors of 2 x 2 x 5 x 59
greySize = 80; %length of grey screen in frames
frameOffset = 0; %number of additional frames before considering response(4 frames = 200ms for Ca kinetics);
movieSize = 1180; %Length of movie in frames
nbins = movieSize/binSize; %make sure this is an integer
chosenMovie = 2; %1=Lab Movie, 2=Nat Movie , if 3 merge 1 and 2
nMovies = 1; % can be 1 or 2
%% Extract and save and load  data
%SpikesSignalFile obtained from running Ca2+ through
                %MLSPike (Deneux, 2016) with the following paramaters:
                %decayTime: 1.2                  # estimated decay time for the indicator utilized (in seconds = 1.2 for gcamp6s)
                %dfOneSpike: 0.2                 # estimated size of df for one single spike (units deltaF/F = 0.2 for gcamp6s)
                %nonLinearities: [0.73, -0.05]   # Amount of indicator non-linearities. For GCaMP6s: [0.73, -0.05]
                %sigmaTune: []                   # A priori level of noise. If left empy, MLspike estimates it from the data
                %f0Drift: [0.00000000000001]     # if this parameter is not set, the algorithm assumes that the baseline remains flat; use 0.00000000000001
%saved in "data"
%data contains inferred spikes in structure of rois x trials x bins, where bin is 1s time bins of the movie

%% calculate whole movie decoding FOR CHOSEN MOVIE
sampleN = min(sessResponderN);
iterations =  100;%iteration_paramater;
[decoderAccuarcy] = MLEDecoderNatStimLOO(data,sampleN, iterations);




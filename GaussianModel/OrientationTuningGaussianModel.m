%% OrientationTuningGaussianModel
% assumes a Gaussian noise model and outputs a p(spike)-based orientation
% tuning curve given subthreshold tuning and trial-to-trial variability for
% CTR and FR groups
%% Load variables
close all;
clear all;

% mean and std (trial-to-trial variability) curves  from data on orientation-evoked subthreshold depolarization relative to rest (group (con, fd) x orientation)
% values for pairs of angles equidistant from preferred angles have been
% averaged to impose symmetry
meanSubThreshTuningCurve = [20.8269518100000,20.9097085650000,21.6898441500000,22.8925377800000,21.6898441500000,20.9097085650000,20.8269518100000;17.6162951025045,17.6862941801951,18.3461650441938,19.3634529362131,18.3461650441938,17.6862941801951,17.6162951025045];
stdSubThreshTuningCurve = [4.22710286480114,4.22710286480114,4.22710286480114,4.22710286480114,4.22710286480114,4.22710286480114,4.22710286480114;4.52360577107142,4.52360577107142,4.52360577107142,4.52360577107142,4.52360577107142,4.52360577107142,4.52360577107142];
% meanEffDistToSpike from data [con, fd]
SpikeThreshold = [32; 27]; %relative to resting potential 

%% calculate pspike
uSpikeThreshold= SpikeThreshold.*ones(2,nOrientations);
uResponse = meanSubThreshTuningCurve;
varResponse = stdSubThreshTuningCurve.^2;
sigmaResponse = sqrt(varResponse);

Zscore = (uSpikeThreshold-uResponse)./sigmaResponse;    
pSpike = normpdf(Zscore); 
pSpike = pSpike./max(pSpike,[],2); %normalize to max


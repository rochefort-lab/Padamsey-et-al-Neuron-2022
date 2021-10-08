function [decoderAccuarcy] = MLEDecoderNatStimLOO(data,sampleN, iterations)
%MLE decoder for responses to natural stimuli using leave one out method
% data contains inferred spikes in structure of rois x trials x bins, where bin is 1s time bins of the movie
%data structure was were obtained from running Ca2+ through MLSPike (Deneux, 2016) which has been integrated with the Rochefort Lab pipeline with the following paramaters:
                %decayTime: 1.2                  # estimated decay time for the indicator utilized (in seconds = 1.2 for gcamp6s)
                %dfOneSpike: 0.2                 # estimated size of df for one single spike (units deltaF/F = 0.2 for gcamp6s)
                %nonLinearities: [0.73, -0.05]   # Amount of indicator non-linearities. For GCaMP6s: [0.73, -0.05]
                %sigmaTune: []                   # A priori level of noise. If left empy, MLspike estimates it from the data
                %f0Drift: [0.00000000000001]     # if this parameter is not set, the algorithm assumes that the baseline remains flat; use 0.00000000000001
%sampleN is how many rois should be used for decoding
%iterations is how many iterations should be run; decoding is averaged
%across iterations

X = data; 
[totalrois, ntrials, nbins] = size(X); % a bin is a scene (1s duration of movie)
Xmean = nan(totalrois, ntrials, nbins);
Xstd = nan(totalrois, ntrials,nbins);
proportionCorrect_it = nan(iterations,1); %holds proportion correct responses per iteration (used for bootstrapping) and per spatial frequency

%generate LOO means and stdevs (ie leave target trial out);
for itrial = 1:ntrials
    idx = [1:ntrials];
    idx(idx==itrial) = [];
Xmean(:,itrial,:) = squeeze(nanmean(X(:,idx,:),2)); %roi x ntrials x bin
Xstd(:,itrial,:) = squeeze(nanstd(X(:,idx,:),[],2)); %roi x ntrials x bin
end

%remove 0 from STD and replace with NaN; 0 std means all trial values are
%0, impossible to calculate probabilities
Xstd(Xstd==0)=NaN;

%X = permute(X,[1,3,2]);% rearrange so nrois x bins x nreps
p_roi = nan(sampleN,nbins,nbins,ntrials);

for it = 1:iterations
    %subsample
    roiSet = datasample([1:totalrois]',sampleN,1,'Replace',false);
    X_sample = X(roiSet,:,:);
    Xmean_sample = Xmean(roiSet,:,:);
    Xstd_sample = Xstd(roiSet,:,:);
    nCorrect = 0; %zero count
    nTotal = 0; %zero count

    for iroi = 1:sampleN %do 1 roi at a time
        
        x  = squeeze(X_sample(iroi,:,:));  %tile values such that a given row a has all bins for a trial, columns are identical
        x = repmat(x,1,1,nbins);
        x = permute(x,[3,2,1]); %nbins x nbins x trials; rows are identical; col holds all bin values for a trial 
        
        mu = Xmean_sample(iroi,:,:); %tile means such that a given row a has identical means for bin a
        mu = repmat(mu,1,1,1,nbins);
        mu = squeeze(mu);
        mu = permute(mu,[2,3,1]);  %nbins x nbins x trials; columns are identical; rows hold all bin means for a given trial
        
        sigma = Xstd_sample(iroi,:,:); %tile stds such that a given row a has identical stds for bin a
        sigma = repmat(sigma,1,1,1,nbins);
        sigma = squeeze(sigma);
        sigma = permute(sigma,[2,3,1]);  %nbins x nbins x trials; columns are identical; rows hold all bin means for a given trial
        
        
        
        Zscore = (x-mu)./sigma;    
        p = normpdf(Zscore);    

    end
    p_roi(iroi,:,:,:) = p;
    likelihood = squeeze(prod(p_roi,1,'omitnan')); %multiply over rois
    
    MLE = (likelihood== max(likelihood,[],1)); %find MLE; nbins x nreps
    MLE = MLE./nansum(MLE,1); %incase of multiple maxes, standardizes.
    totalMax = sum(MLE);

    
    diaganol = diag(MLE); %correct hits
    nCorrect = nCorrect + nansum(diaganol./totalMax');
    nTotal = nTotal+ nbins;
    proportionCorrect_it(it) = nCorrect./nTotal; %collapse across orientaitons; ans = 1x nspatfreqs
 
end
    

decoderAccuarcy = nanmean(proportionCorrect_it);

end
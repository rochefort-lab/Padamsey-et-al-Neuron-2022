function [decoderAccuarcy] = MLEDecoderOrientationLOO(data,N,iterations)
%MLEDecoder for Orientations using leave one out method


%data is a rois x orients x spatfreqs x trials matrix that holds
%single-value stimulus responses across trials. Matrix has NaN values due
%to unequal trial sizes that are dealt with 

%N is how many rois should be used for decoding
%iterations is how many iterations should be run; decoding is averaged
%across iterations

[nrois, norients, nspatfreqs, ~] = size(data);
proportionCorrect_it = nan(iterations,1); %holds proportion correct responses per iteration (used for bootstrapping) and per spatial frequency

for it = 1:iterations
    roiset = datasample([1:nrois]',N,1,'Replace',true);
    
    if N==1
       roiset = [roiset,roiset]; %to perserve dimensions
    end
    
    roiset_data = data(roiset,:,:,:);
    nCorrect = 0; %zero count
    nTotal = 0; %zero count
    for ispatfreq = 1:nspatfreqs %averages across spatial frequency
        
        Xsp_data = squeeze(roiset_data(:,:,ispatfreq,:)); %holds nrois x norientsxnspatfreq vector of all responses for a given rep
        %find relevant trials (those without all NaNs)
        x = squeeze(Xsp_data(:,1,:)); %holds nrois x nreps vector of all responses for a given rep
        relevantTrials = find(sum(isnan(x))~=N);
        Xsp_data = Xsp_data(:,:,relevantTrials);
        nRelevantTrials = numel(relevantTrials);
                
        for iTrial = 1:nRelevantTrials
            X = squeeze(Xsp_data(:,:,iTrial)); %holds nrois x norients vector for a given rep
            omittedOrientations = find(sum(isnan(X),1)==size(X,1)); %all omissions for when >5% of rois show outlier response
            Xsp_LOO = Xsp_data;
            Xsp_LOO(:,:,iTrial) = NaN; %remove test trials from mean/std
           
            
            Xsp_mean = nanmean(Xsp_LOO,3);
            Xsp_std = nanstd(Xsp_LOO,[],3);
            
            X = repmat(X,1,1,norients); %nrois x norients, x copied  norient times in 3rd dimension
            X = permute(X,[1,3,2]); %change to nrois x norients x nreps; reps occur across columns (ie within a row)
            mu = repmat(Xsp_mean,1,1,norients);   %nrois x norients, x copied  norient times in 3rd dimension
            sigma = repmat(Xsp_std,1,1,norients);  %nrois x norients, x copied  norient times in 3rd dimension
            
            Zscore = (X-mu)./sigma;    
            p = normpdf(Zscore); 
            likelihood = squeeze(prod(p,1,'omitnan')); %norients nreps
            MLE = double(likelihood==max(likelihood,[],1)); %find MLE; norients x nreps
            
            MLE(:,omittedOrientations) = NaN; %ignore omitted orientations
            MLEwrap = MLE(1:norients/2,:)+MLE(1+norients/2:end,:);
            MLEwrap(MLEwrap>0)=1; %force multiple hits along wrapped orientations to 1
            
            diaganol = [diag(MLEwrap(:,1:norients/2));diag(MLEwrap(:,1+norients/2:end))];
            totalMax = sum(MLEwrap,1);
            
            nCorrect = nCorrect + nansum(diaganol./totalMax');
            nTotal = nTotal+ norients - numel(omittedOrientations);
            
             
        end
    end
    
    proportionCorrect_it(it) = nCorrect./nTotal; %collapse across orientaitons; ans = 1x nspatfreqs
end


decoderAccuarcy = nanmean(proportionCorrect_it);


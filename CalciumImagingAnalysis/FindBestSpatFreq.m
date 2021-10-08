function     [pref_spatfreq, numResponsiveSF, pref_dir] = FindBestSpatFreq(sig_spatfreqs, mean_dF_response, mean_response_std)
%  finds best spat freq in data set, defined as being the spatial frequency
%  containing the largest max response and being significant

mean_dF = squeeze(mean_dF_response);
maxdF_sig_sf = squeeze(max(mean_dF(:,sig_spatfreqs),[],1)); %local max mean for all sig spat freqs
maxdF_sig_sf(maxdF_sig_sf<0)=0; %exclude any negative responders
numResponsiveSF = sum(maxdF_sig_sf>0); %this calculates number of responsive spatial frequencies. a spatfreq needs to be significant and have a >0 mean response to be considered 

if numResponsiveSF>0
    %dF method
    prefMean = max(maxdF_sig_sf); %max dF response
    pref_spatfreq = sig_spatfreqs(find(maxdF_sig_sf == prefMean,1,'last'));
    pref_dir = find(mean_dF(:,pref_spatfreq)== prefMean,1,'last');
    norients = size(mean_dF,1);
    orth_orient = mod(pref_dir+norients/4,norients);
    if orth_orient == 0
    orth_orient = 1;
    end    
else
    pref_dir = NaN;
    pref_spatfreq = NaN;
end


end



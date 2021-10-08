function [X] =  ProcessCalciumSignalsforGratings(signal, tt, stimtime, stimID) 
%Input:calcium data (signal), time vector (tt), stimulus presentation times 
%(stimtime) and stimulus angle IDs (stimIDs) Relevant outputs are as follows: 

%Output X:
%1 X{1,1} = mat2cell(signalAvgT,nrois,norients,nspatfreqs,avg_framelength*3,3); %contains averaged dF data for each stim condition (signalAvgT)
%2 X{1,2} = mat2cell(deltaBIC_prefspatfreq',1); %contains best deltaBIC value for each ROI 
%3 X{1,3} = mat2cell(nrois,1); %contains total number of ROIs
%4 X{1,4} = mat2cell(pref_spatfreq,nrois); %contains preferred spatial frequencies
%5 X{1,5} = mat2cell(osi,nrois); %contains osi
%6 X{1,6} = mat2cell(pref_orient,nrois); %contains pref orient
%7 X{1,7} = mat2cell(orientations,1); %contains orientations of stim
%8 X{1,8} = mat2cell(spatfreqs,1); %contains spatfreqs of stim
%9 X{1,9} = mat2cell(mean_response,nrois,norients,nspatfreqs); %contains mean signal responses (mean, std, n)
%10 X{1,10} = mat2cell(RTSCfit_orientation,nrois,4); holds values of fit,3rd value = tuning width
%11 X{1,11} = mat2cell(trimmed_response_across_trials,nrois,norients,nspatfreqs,ntrials); %holds trimmed responses acorss trials (trimmed mean specified by sigma)


%Input Format:
%signal = [rois x time x trial) 
%spikes = [rois x time x trial)
%tt = [1 x time x trial]
%stimID = [trial x 1] cell
%stimtime = [trial x 1]

predefinedsf = [0.02, 0.04, 0.08, 0.16, 0.32]; %if undefined finds sf from data
trim_sigma = 2.5; %to assess signficicance of average response
deltaBICthreshold = 10; %model significance for responsivity is based on BIC differences of >= 10 
% Base definitions

imagingfreq = 1/tt(2);
[nrois, ~, ntrials] = size(signal);

% Find values for unique orientations and spatial frequencies used in the experiment
uniquestim = unique([stimID{:}]);

for ngratings=1:numel(strfind([uniquestim{:}],'F')) %gratings defined by presesnce of F in uniqueStim (and therefore stimid)
    o(ngratings) = str2num(uniquestim{ngratings}(3:5)); %position 3:5 yield orientation angle
    sf(ngratings) = str2num(uniquestim{ngratings}(21:24)); %position 21:24 yield spatial frequency
end

orientations = unique(o); %holds a list of unique orientations used in ascending order
norients = numel(orientations);

if ~isempty(predefinedsf)
    spatfreqs = predefinedsf; % hack to ensure equal sized matricies across exps where different numbers of spatfreqs were used
else
    spatfreqs = unique(sf); %holds a list of unique spatial frequencies used in ascending order
end
nspatfreqs = numel(spatfreqs);


% Parse - creates parsed_data, a 5D cell matrix (roi,orientation,spatfreq, trial, data) of signals during stim presentation
parsed_data_stim = zeros(norients,2,ntrials); %initialize parsed_data stim (will hold stimid across trials) 2nd dimension toggles between orientation and spatial freq information

for itrial = 1:ntrials
    %Get trial stim info
    trial_stims = zeros(numel(stimID{itrial}),2); %initialize a matrix that holds x(1)=vis stim (orienation,grey,black), x(2)=spatial freq for itrial
    
    for istim = 1:numel(stimID{itrial}) %fills trial_stims;
        if contains(stimID{itrial}{istim},'F')%if stim is a grating
            trial_stims(istim,1) = str2num(stimID{itrial}{istim}(3:5));
            trial_stims(istim,2) = str2num(stimID{itrial}{istim}(21:24));
        end
    end
    
    parsed_data_stim(:,:,itrial) = trial_stims(trial_stims(:,2)>0,:); %holds trial_stims data across trials, but removes zeros
    trialspatfreq = parsed_data_stim(1,2,itrial); %holds spatial frequency for trial
    ispatfreq = find(spatfreqs==trialspatfreq); %holds spatial frequency posistion for addressing arrays
    [~, order_of_stim] = sort(parsed_data_stim(:,1,itrial));%re-orders orienations in parsed_data_stim in ascending order
    %order_of_stim holds index of positions of original orientations
    
    stim_indicies = find(trial_stims(:,2)); %find position of stimuli presented, so can retreive start and end frames for each stimulus: see next 2 lines
    stim_startframe = ceil(stimtime{itrial}(stim_indicies)*imagingfreq); %finds starting frame for each stim presentation
    stim_endframe = floor(stimtime{itrial}(stim_indicies+1)*imagingfreq); %finds ending frame for each stim presentation
   %driftDuration * imagingfreq; % calculates average duration of stim presentation
    %initialize parsed_data array now that avg_framelength is known
    if itrial == 1
        avg_framelength = floor(mean(stim_endframe - stim_startframe)); 
        parsed_data = NaN(nrois,norients,nspatfreqs,ntrials,avg_framelength*3); %initialize parsed_data, holds dF; %initialize parsed_data, holds dF

    end
    

    %for each roi, extracts relevant portion of recording, and orders it in
    %parsed_data matrix(orientation,trial,spatfreq,roi) in accordance to
    %orienation and spatfreq
    
    for iroi=1:nrois
        %zero  the trace, assume smallest 5% of response is baseline
        norm_signal = signal(iroi,:,itrial)-quantile(signal(iroi,:,itrial),0.05);
        for iorient=1:norients %for each unique stim, pulls relevant raw data from the matrix signals for the duration of stimulus presentation + and - the avg_num of frames (to get baseline, and decay)
                s = norm_signal(stim_startframe(order_of_stim(iorient))-avg_framelength:stim_startframe(order_of_stim(iorient))+2*avg_framelength-1);
                if ~isempty(ispatfreq)%might not work if trial spatfreq is not on predefined spatfreqs listed 
                    parsed_data(iroi,iorient,ispatfreq,itrial,:)= s; %at2cell(s,1,numel(s));
                end
        end
    end
end
allTrialSpatFreqs = squeeze(parsed_data_stim(1,2,:)); %gives ntrials x 1 array of spatfreq in trial



%for  rois, and for a given orientation and spat freq, avg all
%trials and store avg, std, n
%initialize matricies
signalAvgT = zeros(nrois,norients,nspatfreqs,avg_framelength*3,3); %initialize X, holds mean, std, n (toggeled by last dimension (1,2,3))
mean_response = zeros(nrois,norients,nspatfreqs);
mean_response_std = zeros(nrois,norients,nspatfreqs); %holds mean/SE response of stim and poststim period, minus, baseline
mean_response_n = zeros(nrois,norients,nspatfreqs); %holds mean/SE response of stim and poststim period, minus, baseline


%initialize matricies for pref_spatfreq, and for osi and dsi analysis
pref_spatfreq = nan(nrois,1);
numResponsiveSF = zeros(nrois,1);

osi = nan(nrois,1);
pref_orient = nan(nrois,1);
prefDirInt = nan(nrois,1);
deltaBIC=nan(nrois,nspatfreqs);%holds BIC per spatfreq
deltaBIC_prefspatfreq = NaN(nrois,1); %holds final deltaBIC value for pref spat freq


response_across_trials = nan(nrois,norients,nspatfreqs,ntrials);
baseline_across_trials = nan(nrois,norients,nspatfreqs,ntrials);
trimmed_response_across_trials = nan(nrois,norients,nspatfreqs,ntrials);

RTSCfit_orientation = nan(nrois,4); %holds Rp, theta, sigma, C for gaussian fit

%loop for calculating signalAvg (dF vs t) and stimulus response (dF)
%averages

%for calculating baseline before calcium imaging response
windowSize = 40; %in frames (2s)
usedWindowSize = 40; %frames from end of baseline window
for iroi = 1:nrois
    
    for ispatfreq=1:nspatfreqs
        relevantTrials = (allTrialSpatFreqs==spatfreqs(ispatfreq)); %relevant trials are those with given spatial freq
        nrelevantTrials = sum(relevantTrials);
        sq = squeeze(parsed_data(iroi,:,ispatfreq,relevantTrials,:)); %holds all trials for a given stim condition
        
        %find max_response_period based on global avg response
        global_avg = squeeze(nanmean(nanmean(sq,1),2));
        smoothed_global_avg =  movmean(global_avg(avg_framelength+1:end),avg_framelength,'Endpoints','discard'); %done for the purposes of finding best avg_framelength response period
        [~,max_smoothed_global_avg] = max(smoothed_global_avg); %finds best avg_framelength response period
        max_response_period_dF = max_smoothed_global_avg+avg_framelength;
        
        for iorient = 1:norients
            sq_orient = squeeze(sq(iorient,:,:)); %holds all trials for agiven stim condition
            n_orient = size(sq_orient,1) - sum(isnan(sq_orient(:,1)));
              
                %avg, std, n across time 
                signalAvgT(iroi,iorient,ispatfreq,:,1) = nanmean(sq_orient); %holds mean vs t
                signalAvgT(iroi,iorient,ispatfreq,:,2) = nanstd(sq_orient); %holds std vs t
                signalAvgT(iroi,iorient,ispatfreq,:,3) = n_orient*ones(1,avg_framelength*3); %holds N vs t
                                       
                baselineWindows =  movmean(sq_orient(:,avg_framelength-windowSize+1:max_response_period_dF),windowSize,2, 'Endpoints', 'discard');
                [~, loc_min_baseline] = min(baselineWindows,[],2); %local reference
                loc_min_baseline_start = loc_min_baseline+avg_framelength-windowSize; %change reference
                loc_min_baseline_start = loc_min_baseline_start+(windowSize-usedWindowSize);
               
                %take windowSized baselinbe from calcualted location
                b = nanmean(sq_orient(:,loc_min_baseline_start:loc_min_baseline_start+usedWindowSize-1),2);
               
                baseline_across_trials(iroi,iorient,ispatfreq,1:nrelevantTrials) = b;
                sq_response_across_trials =  mean(sq_orient(:,max_response_period_dF:max_response_period_dF+avg_framelength-1),2)-squeeze(baseline_across_trials(iroi,iorient,ispatfreq,1:nrelevantTrials));
                                
                response_across_trials(iroi,iorient,ispatfreq,1:nrelevantTrials) = sq_response_across_trials;                
                trimmed_responses = trim(sq_response_across_trials,1,trim_sigma); %trim outlier responses, defined by >2.5x stdev of response variability 
                trimmed_response_across_trials(iroi,iorient,ispatfreq,1:nrelevantTrials) = trimmed_responses;                
                
                %TRIMMED means, std, n
                mean_response(iroi,iorient,ispatfreq) = nanmean(trimmed_responses); %nanmean(sq_response_across_trials);
                mean_response_n(iroi,iorient,ispatfreq) = sum(~isnan(trimmed_responses));%n_orient; %holds n associated with mean
                mean_response_std(iroi,iorient,ispatfreq) = nanstd(trimmed_responses); %nanstd(sq_response_across_trials); %holds std associated with mean
            end
        
        %test significance based on gaussian fits
           if iroi>2 %roi 1 and 2 are background and npil
           [deltaBIC(iroi,ispatfreq)] = GaussianFit_Significance(squeeze(trimmed_response_across_trials(iroi,:,ispatfreq,1:max(mean_response_n(iroi,:,ispatfreq))))');
           end
    end

    sig_spatfreqs = find(deltaBIC(iroi,:)>=deltaBICthreshold); %holds all sig spatfreq (from stats test)
    if ~isempty(sig_spatfreqs) && iroi>2 %if correct_npil is true and roi correlation with npil exceeds defined threshold, then omit roi
        [pref_spatfreq(iroi), numResponsiveSF(iroi), prefDirInt(iroi)] = FindBestSpatFreq(sig_spatfreqs, mean_response(iroi,:,:), mean_response_std(iroi,:,:)); %find which spat freq best evokes max sig activity in neuron and use it for further analysis
        if numResponsiveSF(iroi)>0
            [pref_orient(iroi),osi(iroi)] = getOSI(mean_response(iroi,:,pref_spatfreq(iroi)),orientations);
            RTSCfit_orientation(iroi,:) = GaussianFit_Orientation(iroi,squeeze(trimmed_response_across_trials(iroi,:,pref_spatfreq(iroi),1:max(mean_response_n(iroi,:,pref_spatfreq(iroi)))))');
            deltaBIC_prefspatfreq(iroi) = deltaBIC(iroi,pref_spatfreq(iroi));
        end
    end
    
end

    
%OUTPUT 
X{1,1} = mat2cell(signalAvgT,nrois,norients,nspatfreqs,avg_framelength*3,3); %contains averaged dF data for each stim condition (signalAvgT)
X{1,2} = mat2cell(deltaBIC_prefspatfreq',1); %contains best deltaBIC value for each ROI 
X{1,3} = mat2cell(nrois,1); %contains total number of ROIs
X{1,4} = mat2cell(pref_spatfreq,nrois); %contains preferred spatial frequencies
X{1,5} = mat2cell(osi,nrois); %contains osi
X{1,6} = mat2cell(pref_orient,nrois); %contains pref orient
X{1,7} = mat2cell(orientations,1); %contains orientations of stim
X{1,8} = mat2cell(spatfreqs,1); %contains spatfreqs of stim
X{1,9} = mat2cell(mean_response,nrois,norients,nspatfreqs); %contains mean signal responses (mean, std, n)
X{1,10} = mat2cell(RTSCfit_orientation,nrois,4); %holds values of fit,3rd value = tuning width
X{1,11} = mat2cell(trimmed_response_across_trials,nrois,norients,nspatfreqs,ntrials); %holds trimmed responses acorss trials (trimmed mean specified by sigma)



end
%% CalciumImagingAnalysisforGratings
% analyzes calcium imaging data for grating-evoked responses, mainly
% orientation tuning responses (OSI, OSIGaussian), but also orientation
% decoding


%% Run
clear all
close all;
% first run doProjFuncAna(data_dir, projName, 'deltafspikes', 'Zahid_FDepAnalysisSpikes', {'drift'}, [], false, 20, true, [], true, 'signalOnly');
% Inputs: doProjFuncAna(data_dir, projName, signalName, funcStr, groupNames, ...
%    collapseFunc, flagCellMode, analysis_fs, includeVSM, batchTiffStacks,...
%    includeBG, auxSignal) %auxSignal = 'signalOnly' bypasses wheel=
%this analysis pushes data through the Rochefort pipeline, finally calling on ProcessCalciumSignalsforGratings
X = ProcessCalciumSignalsforGratings; % outputs X
%DATA structure returned from ProcessCalciumSignals:
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
 
%% Extract variables from main data

%re-define stimuli for ease
signal = cell2mat([data{:,1}]');
[totalrois, norients, nspatfreqs, trial_duration, ~] = size(signal);
signal_means = cell2mat([data{:,9}]');
OrientationRTSC = cell2mat([data{:,10}]');
OSIGaussian = OrientationRTSC(:,3);
deltaBIC_prefspatfreq = cell2mat([data{:,2}])';
prefSpatFreq = cell2mat([data{:,4}]');
OSI = cell2mat([data{:,5}]');
prefOrient = cell2mat([data{:,6}]');
nrois = cell2mat([data{:,3}]');
orientations = unique(cell2mat([data{:,7}]'));
%ensure it is a row vector
if size(orientations,1) ~=1
    orientations = orientations';
end
spatfreqs = unique(cell2mat([data{:,8}]'));
uid = sesname; 
nsessions = length(uid);
uid = string(uid);
%% Calculate tuning curve (per roi)
tuningCurve = nan(totalrois,norients/2); %single roi tuning cover with normalized mean
assessedRois = [find(~isnan(prefSpatFreq))]';
for iassessedRoi = 1:numel(assessedRois)
    iroi = assessedRois(iassessedRoi);
    %Tuning curve
    u = mean([signal_means(iroi,1:norients/2,prefSpatFreq(iroi));signal_means(iroi,norients/2+1:end,prefSpatFreq(iroi))],1); %FOLD DIRECTIONS
    u = u./max(u); %normalize to max
    tuningCurve(iroi,:)=CenterTuningCurve(u);%align to center

end
%% Assemble AN_roi_raw matrix, holds variables of interest per roi
AN_roi_raw = [OSI, OSIGaussian];
AN_labels = strings(1,size(AN_roi_raw,2));
AN_labels(1) = 'OSI'; AN_labels(2) = 'OSIGaussian';
%% Create experimental tags
%get experimental IDs for each roi
nuids = length(sesname);
condition_tag = strings(totalrois,1); %holds roi condition (ie fd, fdr)
uid_condition = strings(nuids,1); %gets condition IDs, uses to fill condition_tag
%get uid_condition from yaml incase of any changes
for iuid = 1:nuids
metaData_file = strcat(fullfile('C:\Users\Zahid\Documents\GitHub\NRLabMetadata\sessions\Zahid',sesname(iuid)),'.yaml'); %holds session metadata
sessInfo = ReadYaml(metaData_file{1});
uid_condition(iuid) = sessInfo.recordings.(animal{iuid}).calciumCells{1};
end

roi_tag = zeros(totalrois,1);  %holds roiIDs for each exp
animal_tag = strings(totalrois,1);
uid_animal = string(animal); %animal is the original variable, will be replaced with uid_animal
[animals ~] = unique(uid_animal,'stable');
nanimals = numel(animals);
uid_tag = strings(totalrois,1);
ses_tag = strings(totalrois,1);
conditions = unique(uid_condition,'stable');
animal_condition = uid_condition(animal_loc);
nconds = length(conditions);

for i=1:nsessions%size yields number of experiments
    start = min(find(condition_tag=="")); %finds first empty spot in condition_tag
    condition_tag(start:start+nrois(i)-1) = uid_condition(i); %fill tag with exp id found in markedtype vector
    %     chronicroi_tag(start:start+nrois(i)-1) = chronicroi_expID(i); %fill tag with chronic roi id
    animal_tag(start:start+nrois(i)-1) = uid_animal(i); %fill tag with animal id found in markedtype vector
    uid_tag(start:start+nrois(i)-1) = uid(i); %fill tag with animal id found in markedtype vector
    ses_tag(start:start+nrois(i)-1) = sesname(i); %fill tag with animal id found in markedtype vector
    roi_tag(start:start+nrois(i)-1) = [1:nrois(i)];
end
uid_conditionsInt = zeros(nuids,1);
uid_tag_Int = zeros(totalrois,1);
for iuid = 1:nuids    
uid_conditionsInt(iuid) = find((conditions==uid_condition(iuid)));
uid_tag_Int(uid_tag==uid(iuid)) = iuid;
end
%% Only use significant rois (BIC>=10)
deltaBIC_threshold = 10;
sigrois = deltaBIC_prefspatfreq<deltaBIC_threshold;
roi_filter_bin = useable_roi_bin & sigrois; %roi 1 and 2 are background and npil, true rois start after 
roi_filter = double(roi_filter_bin); %converts to double for purposes of replacing 0s with NaNs
roi_filter(roi_filter==0) = NaN; %replace 0s with NaNs
AN_roi = AN_roi_raw./roi_filter;
%% Average within rois within a session (uid) and condition (conds)
%initalize variables
AN_session = zeros(nuids,size(AN_roi,2));
AN_condition = zeros(nconds,size(AN_roi,2));

%include response fraction
if ~sum(contains(AN_labels,'responseFrac'))
    AN_labels(end+1) = 'responseFrac';
end
AN_session(:,end+1)=zeros(nuids,1);
AN_condition(:,end+1)=zeros(nconds,1);

%include nrois
if ~sum(contains(AN_labels,'nSelectedRois'))
    AN_labels(end+1) = 'nSelectedRois';
end
AN_session(:,end+1)=zeros(nuids,1);
AN_condition(:,end+1)=zeros(nconds,1);

%for tuning curve calculation per condition
tuningCurve_condition = NaN(nconds,norients/2,4,3); %1st dim: nconds; 2nd dim: orients 3rd dim: mean, SE, N
tuningCurve_session = NaN(nsessions,norients/2); %1st dim: nsessions; 2nd dim: orients

%calculate session averages
nOrigLabels = size(AN_roi,2);
nselectedRois = zeros(nuids,1);
for iuid=1:nuids
    x = uid_tag==uid(iuid);
    AN_session(iuid,1:nOrigLabels) = nanmean(AN_roi(x,:),1);
    AN_session(iuid,nOrigLabels+1) = nansum(roi_filter_bin(x))/(nansum(useable_roi_bin(x))); %nansum(roi_filter_bin(x))/(nansum(usable_roi_bin(x)));%response fraction;
    nselectedRois(iuid) = sum(~isnan(roi_filter(x)));
    AN_session(iuid,nOrigLabels+2) = nselectedRois(iuid);%total rois
    tuningCurve_session(iuid,:) = nanmean(tuningCurve(roi_filter_bin&x,:),1);
end

%calculate condition averages
for icond=1:nconds
    x = condition_tag==conditions(icond); %tags for Rois
    y = uid_condition==conditions(icond);
    AN_condition(icond,:) = nanmean(AN_session(y,:),1);
    %tuningCurve averaged across session
    tuningCurve_condition(icond,:,1) = nanmean(tuningCurve(roi_filter_bin&x,:),1); %holds mean
    tuningCurve_condition(icond,:,2) =  nanstd(tuningCurve_session(y,:),[],1)./sqrt(sum(y)); %holds se
    tuningCurve_condition(icond,:,3) = sum(y); %holds n
end
%% MLE Decoding - Orienation (fixed n; LOO)
decoderAccuarcy_session = nan(nuids,1);
distanceOfError_session = nan(nuids,1);
N = 50; %number of cells to use for a given decoding run
iterations = 100;
%only do first session for animal
[~,firstSession] = unique(uid_animal);
for iuid = firstSession'; %1:nuids
    sprintf('Doing session %d of %d',iuid,nuids)
    trialResponses =  cell2mat(data{iuid,11}'); 
    xall = find(uid_tag==uid(iuid) & useable_roi_bin); %roi_tag index for session + omissions
    x = xall-find(uid_tag==uid(iuid) & roi_tag==1)+1; %subtraction to normalize roi indicies to the uid session
    [decoderAccuarcy_session(iuid)] = MLEDecoderOrientationLOO(trialResponses(x,:,:,:,:),N,iterations);
end

%calc condition avg
decoderAccuarcy_condition = nan(nconds,1);
distanceOfError_condition = nan(nconds,1);
for icond=1:nconds
    y = uid_condition==conditions(icond); %tags for sessions
    decoderAccuarcy_condition(icond) = nanmean(decoderAccuarcy_session(y)); %weight by nrois
end

function [X,samplingFreq,abfTimeStamp] = getABFfiles(RecordingFilePaths,DesiredSamplingFreq)
%load ABF files, primary channel is extracted and stored in X along with
%timestamp of recording. Recordings are downsampled to DesiredSamplingFreq;
%actual sampling frequency is stored in samplingFreq
X = [];
nFiles = length(RecordingFilePaths);
abfTimeStamp = nan(nFiles,1);
iCounter = 1;
        %intrinsicProperties_recording = nan(nFiles,nProperties);
        for iFile = 1:nFiles %get data

            DataStruct = import_wcp(char(RecordingFilePaths(iFile)));
            downSampleFactor = round(1/DataStruct.t_interval/DesiredSamplingFreq);
            samplingFreq = 1/DataStruct.t_interval/downSampleFactor;
            D = downsample(DataStruct.S{1},downSampleFactor);%only take primary channel
         
            %get time stamp
            t = DataStruct.time(strfind(DataStruct.time,' ')+1:end)';
            t_int=str2double(string(t)); %detect numbers
            t_int(isnan(t_int))=[];
            
            %store
            for iTrace = 1:nTraces
            X(:,iCounter) = D(:,iTrace);
            abfTimeStamp(iCounter) = (t_int(1)*10+t_int(2))*3600+(t_int(3)*10+t_int(4))*60+(t_int(5)*10+t_int(6));
            iCounter = iCounter+1;            
            end
        end
        
        
end


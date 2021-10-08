function [rearranged , index] = CenterTuningCurve( data )
%centers tuning curves with respect to their preferred orientation
%data is trial x orienation matrix (works with any number of orientations)
%rearranged = data is centered with preferred response in middle; for 6 orientations data is centered as -60 -30 0 30 60 90 degrees from preferred; where 0 is preferred
%index = an array of values such that data(index) = rearranged

 rearranged = NaN(size(data));
 u =  mean(data,1);
 norients = numel(u);
 disp = mod(round(norients/2)-find(u==max(u),1,'first'),norients);%displacement from center
 index = nan(1,norients);
    for i = 1:norients
        index(i) = mod(i+disp,norients);
        if index(i) ==0
        index(i) = norients;
        end
     rearranged(:,index(i)) = data(:,i);
    end


end


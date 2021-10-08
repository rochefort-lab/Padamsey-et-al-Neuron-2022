function [pref_orient,OSI] = getOSI(data,orientations)
%calculates OSI as 1-circular variance
%data is a row vector with mean values for direction, orientaitons is a row vector that holds
%all directions (not yet collapsed into orientations)
if size(data,2)==numel(orientations)
data = mean(reshape(data,numel(data)/2,2)');%average same orientations;

end
data(data<0)=0; %zero any negative values so that max(osi) == 1

%calculate angle in radians
angle_rad = degtorad(orientations(1:length(orientations)/2));
    % Calculate osi by circular variance method
     vector = sum(data.*exp(2*1i*angle_rad))/sum(data);
     pref_orient = radtodeg(mod(angle(vector), 2*pi))/2;% Move from (-pi,pi] to [0,2pi) representation, divide by 2 as only require half the space for orientation
     OSI = abs(vector);
%      
%      if OSI>1
%          OSI=1
%      end    
        
end
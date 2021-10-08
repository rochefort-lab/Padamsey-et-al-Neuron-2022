function [data] = trim(data,dim_to_trim,sigma,centralTendency, groupIdx)
%function for trimming outlier responses

% Data is 2D
%dim_to_trim = 1; trim along col
%dim_to_trim = 2; trim along row
%groupIdx denotes subgroups of data
if nargin<3
    sigma = 2.5;
end
nn
if nargin<4
    centralTendency = 'mean';
end

if nargin<5
    groupIdx = ones(size(data,dim_to_trim),1);
end

if size(groupIdx,dim_to_trim)==1
    groupIdx = groupIdx';
end


nGroups = numel(unique(groupIdx));

for iGroup = 1:nGroups
           targetIdx = (groupIdx==iGroup);
           
           if dim_to_trim==1
           y = data(targetIdx,:);
           else
           y = data(:,targetIdx);
           end
          if contains(centralTendency,'mean')
           m = nanmean(y,dim_to_trim);
           err = sigma*nanstd(y,[],dim_to_trim);
          else
           m = nanmedian(y,dim_to_trim);
           stdev = nanmedian(abs(y-m),dim_to_trim);
           err = sigma*stdev;
          end
          
           pos_err = m+err;
           neg_err = m-err;
           data(data>pos_err & targetIdx )=NaN; %remove max element
           data(data<neg_err & targetIdx)=NaN; %remove min element
          
end
end


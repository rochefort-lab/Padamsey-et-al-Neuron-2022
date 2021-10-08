function [ X ] = roundToNearest(X,m,direction)
%round array X to the nearest m

if nargin<3
    direciton = 'normal'
end

if contains(direction,'up')
X = round(ceil(X/m))*m;
    
else
X = round(X/m)*m;
end
    
end


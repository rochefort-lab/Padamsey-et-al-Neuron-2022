function [RTSCfit] = GaussianFit_Orientation(iroi,data,clampC,plotData,sigmaInputSet)
%% Fit Gaussian
%X is a matrix of input ntrials x norients for a given roi
%fits gaussian returning values for R (amplitude at pref) T (pref angle) S
%(width of gaussian as defined by 1 std) and C (offset) (returned in
%RTSCfit)

if nargin<3
    clampC = true; %better fit achieved with C is clamped
end

if nargin<4
   plotData=false; 
end

if nargin<5
   sigmaInputSet=[]; %initialization of sigma for gaussian fits can be user-defined
end


[~, norients] = size(data);

%fold data (ie average opposite directions for the same orientation)
X = [data(:,1:norients/2);data(:,norients/2+1:end)];

%step size
alpha_u_theta = 1;
alpha_delta_theta = 0.01;
alpha_sigma = 1;
alpha_C = 0.01;
alpha_Rp = 0.01;

%initialize paramaters
sigmaMin = 10; %min value of sigma (set to 10)

if ~isempty(sigmaInputSet)
sigma_set = [1:8]*(360/(norients*2));  %test 6 initial values of sigma spaced at 1/2 interval of orientations
else
sigma_set=sigmaInputSet;   
end

C_start = 0;

MaxC = 1000; %max value of C if unclamped (default 0)
mean_orientationR = nanmean(X,1);
MaxR = max(mean_orientationR);
Rp_start = MaxR;
theta = linspace(0,180-(180/(norients/2)),norients/2);
u_theta_start =  theta(find(mean_orientationR == max(mean_orientationR),1,'first')); %initialize mean theta as theta with max response


nsigmas = numel(sigma_set);
SSE = zeros(nsigmas,1); %sum of squares
modelParameters = zeros(nsigmas,4); %holds model paramaters for each sigma tested

for isigma = 1:nsigmas
    
    %re-initialize
    u_theta = u_theta_start;
    Rp = Rp_start;
    C = C_start;
   
    for i = 1:1000 %number of iterations
        
        %Gradient descent
        %center theta
        delta_theta = theta-u_theta;
        circ_delta_theta = [delta_theta;delta_theta-180;delta_theta+180]; %wraps around 0-90 for orientation
        delta_theta = abs(min(abs(circ_delta_theta)));
        
        %core functions
        core_gauss = gaussmf(delta_theta,[sigma_set(isigma) 0]);
        core_derivative = X-(Rp*core_gauss+C);
        
        %u_theta
        deriv_delta_theta = theta-(u_theta+alpha_delta_theta);
        circ_deriv_delta_theta = [deriv_delta_theta;deriv_delta_theta-180;deriv_delta_theta+180];
        deriv_delta_theta = abs(min(abs(circ_deriv_delta_theta)));
        deriv_delta_theta = deriv_delta_theta - delta_theta;
        
        d_u_theta =  core_gauss.*sign(deriv_delta_theta).*((delta_theta)./(sigma_set(isigma).^2));
        d_u_theta = nansum(nansum(core_derivative.*d_u_theta))*alpha_u_theta;
        u_theta = u_theta-d_u_theta;
        
        %sigma
        d_sigma =  core_gauss.*(((delta_theta).^2)./(sigma_set(isigma).^3));
        d_sigma = nansum(nansum(core_derivative.*d_sigma))*alpha_sigma;
        sigma_set(isigma) = sigma_set(isigma)+d_sigma;
        
        if sigma_set(isigma)<(sigmaMin)
            sigma_set(isigma)=sigmaMin;
        end    
        
        %C
        if clampC==false
            d_C =  1;
            d_C = nansum(nansum(core_derivative.*d_C))*alpha_C;
            C = C+d_C;

            if C<-MaxC
                C=-MaxC;
            end    

             if C>MaxC
                C = MaxC;
            end    
        end
        %Rp
        d_Rp = core_gauss;
        d_Rp = nansum(nansum(core_derivative.*d_Rp))*alpha_Rp;
        Rp = Rp+d_Rp;
        
        if Rp>MaxR
         Rp = MaxR;
        end
       
        if Rp<MaxR/2
        Rp = MaxR/2;
        end

    end
    model =  Rp*core_gauss+C;
    SSE(isigma) = nansum(nansum((model-X).^2));
    modelParameters(isigma,:) = [Rp,u_theta,sigma_set(isigma),C];

end

%weighting paramaters across models in accordance to r2 value to prevent
%overfitting 
    SST = nansum(nansum((X-nanmean(nanmean(X))).^2));%total error
    SSR = SST-SSE; 
    SSR(SSR<0)= 0; %zero models that do worse than mean
    r2 = (SSR./SST).^2;
    if max(SSR)<=0 %suggests mean is best model, chose model with highest sigma
    [RTSCfit] = modelParameters(end,:);
    else %weight the models
    modelWeight = r2./sum(r2);
    [RTSCfit] = sum((modelParameters).*(modelWeight));
     end
% % % % 
% % Graph
if plotData
 theta_cont = [0:150];

        delta_theta1 = theta_cont-RTSCfit(2);
        circ_delta_theta1 = [delta_theta1;delta_theta1-180;delta_theta1+180]; %wraps around 0-90 for orientation
        delta_theta1 = abs(min(abs(circ_delta_theta1)));
        
    gaussianfit = RTSCfit(1)*gaussmf(delta_theta1,[RTSCfit(3) 0])+RTSCfit(4);

%%scatterX and scatterY shape data matrix to plot indivusal points on the gaussian    
scatterX = ones(size(X)).*theta;
scatterX = reshape(scatterX,1,[]);
scatterY = reshape(X,1,[]);

 figure; plot(theta,mean_orientationR); hold on; plot(theta_cont,gaussianfit); hold on; scatter(scatterX,scatterY);
title(iroi);

end

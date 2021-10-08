function [deltaBIC] = GaussianFit_Significance(X)
%% Fit Gaussian
%X is a matrix of input ntrials x norients for a given roi
%function fits gaussian to orientation responses and returns BIC difference between gaussian fit and a single
%line across zero; a deltaBIC>=10 strongly favours the gaussian fit model and suggests the roi is visually responsive
[~, norients] = size(X);
%fold data
mean_orientationR = nanmean(X,1);
std_orientationR = nanstd(X,[],1);
[MaxR,prefOrient] = max(mean_orientationR);

%step size
alpha = 0.01;
alpha_u_theta = 1;
alpha_delta_theta = alpha;
alpha_sigma = 1;
alpha_C = alpha;
alpha_Rp = alpha;
alpha_Rn = alpha;
%initialize paramaters
%sigma is standard deviation of tuning
%C is offset
%Rp is amplitude of preferred
%Rn is amplitude of orthoganol
%theta is the current angle

sigmaMin = 10; %min value of sigma; set to 10
sigma_set = [1:8]*(360/(norients*2));  %test 6 initial values of sigma spaced at 1/2 interval of orientations
nsigmas = numel(sigma_set);
C_start = 0; 
clampC = 1; %1 means to clampC at C-start value; c is clamped for calcium imaging data as it provides a better fit
Rp_start = MaxR;
Rn_start = MaxR;

theta = linspace(0,360-(360/norients),norients);
u_theta_start =  theta(find(mean_orientationR == max(mean_orientationR),1,'first')); %initialize mean theta as theta with max response

SSE = zeros(nsigmas,1); %sum of squares
modelParameters = zeros(nsigmas,5); %holds model paramaters for each sigma tested

for isigma = 1:nsigmas
    %re-initialize variables
    u_theta = u_theta_start;
    Rp = Rp_start;
    Rn = Rn_start;
    C = C_start;
    

    for i = 1:1000 %number of iterations
        
        %Gradient descent
        %center theta
        delta_theta1 = theta-u_theta;
        circ_delta_theta = [delta_theta1;delta_theta1-360;delta_theta1+360]; %wraps around 0-90 for orientation
        delta_theta1 = abs(min(abs(circ_delta_theta)));
        
        delta_theta2 = theta+180-u_theta;
        circ_delta_theta = [delta_theta2;delta_theta2-360;delta_theta2+360]; %wraps around 0-90 for orientation
        delta_theta2 = abs(min(abs(circ_delta_theta)));
        
        %core functions
        core_gauss1 = gaussmf(delta_theta1,[sigma_set(isigma) 0]);
        core_gauss2 = gaussmf(delta_theta2,[sigma_set(isigma) 0]);
        core_derivative = X-((Rp*core_gauss1) +(Rn*core_gauss2)+C);
        
        %u_theta
        deriv_delta_theta1 = theta-(u_theta+alpha_delta_theta);
        circ_deriv_delta_theta1 = [deriv_delta_theta1;deriv_delta_theta1-360;deriv_delta_theta1+360];
        deriv_delta_theta1 = abs(min(abs(circ_deriv_delta_theta1)));
        deriv_delta_theta1 = deriv_delta_theta1 - delta_theta1;
        
        deriv_delta_theta2 = theta+180-(u_theta+alpha_delta_theta);
        circ_deriv_delta_theta2 = [deriv_delta_theta2;deriv_delta_theta2-360;deriv_delta_theta2+360];
        deriv_delta_theta2 = abs(min(abs(circ_deriv_delta_theta2)));
        deriv_delta_theta2 = deriv_delta_theta2 - delta_theta2;
        
        d_u_theta =  core_gauss1.*sign(deriv_delta_theta1).*((delta_theta1)./(sigma_set(isigma).^2))+(core_gauss2.*sign(deriv_delta_theta2).*((delta_theta2)./(sigma_set(isigma).^2)));
        d_u_theta = nansum(nansum(core_derivative.*d_u_theta))*alpha_u_theta;
        u_theta = u_theta-d_u_theta;
        
        %sigma
        d_sigma =  core_gauss1.*(((delta_theta1).^2)./(sigma_set(isigma).^3)) + core_gauss2.*(((delta_theta2).^2)./(sigma_set(isigma).^3));
        d_sigma = nansum(nansum(core_derivative.*d_sigma))*alpha_sigma;
        sigma_set(isigma) = sigma_set(isigma)+d_sigma;
        
        if sigma_set(isigma)<(sigmaMin)
            sigma_set(isigma)=sigmaMin;
        end    
        
        %C
        if ~clampC
            d_C =  1;
            d_C = nansum(nansum(core_derivative.*d_C))*alpha_C;
            C = C+d_C;
        end
        %Rp
        d_Rp = core_gauss1;
        d_Rp = nansum(nansum(core_derivative.*d_Rp))*alpha_Rp;
        Rp = Rp+d_Rp;
        
        if Rp>MaxR
         Rp = MaxR;
        end
        
        
        %Rn
        d_Rn = core_gauss2;
        d_Rn = nansum(nansum(core_derivative.*d_Rn))*alpha_Rn;
        Rn = Rn+d_Rn;
        
        if Rn>MaxR
         Rn = MaxR;
        end
        
        if Rn<0
         Rn = 0;
        end
       
    end
    model =  Rp*core_gauss1+Rn*core_gauss2+C;
    SSE(isigma) = nansum(nansum((model-X).^2));
    modelParameters(isigma,:) = [Rp,Rn,u_theta,sigma_set(isigma),C]; %holds model paramaters across sigma tested

end


%calculate nullModel - this is just 0
nullC = 0;

%calculate R2 of models
    SST = nansum(nansum((X-nullC).^2));%total error with respect to a null model
    SSR = SST-SSE; %residual error
    SSR(SSR<0)= 0; %zero models that do worse than mean
    r2 = (SSR./SST).^2;
    if max(SSR)<=0 %suggests mean is best model, chose model with highest sigma
     [RNTSCfit] = modelParameters(end,:);
    else %weight the models in accordance with R2
    modelWeight = r2./sum(r2);
    [RNTSCfit] = sum((modelParameters).*(modelWeight));
    end
    
%calculate SSE based on selected model
        %model SSE
        delta_theta1 = theta-RNTSCfit(3);
        circ_delta_theta = [delta_theta1;delta_theta1-360;delta_theta1+360]; %wraps around 0-90 for orientation
        delta_theta1 = abs(min(abs(circ_delta_theta)));
        
        delta_theta2 = theta+180-RNTSCfit(3);
        circ_delta_theta = [delta_theta2;delta_theta2-360;delta_theta2+360]; %wraps around 0-90 for orientation
        delta_theta2 = abs(min(abs(circ_delta_theta)));

        gaussianfit = RNTSCfit(1)*gaussmf(delta_theta1,[RNTSCfit(4) 0])+RNTSCfit(2)*gaussmf(delta_theta2,[RNTSCfit(4) 0])+RNTSCfit(5);
        GaussianSSE = nansum(nansum((gaussianfit-X).^2));

        %null SSE
        NullSSE = nansum(nansum((X-nullC).^2));

        
%Model selection based on BIC

N = nansum(nansum(~isnan(X)));
if ~clampC
    GaussianBIC = N*log(GaussianSSE/N)+5*log(N); %for AIC penalty is 2k vs k*ln(N)
else
    GaussianBIC = N*log(GaussianSSE/N)+4*log(N); %for AIC penalty is 2k vs k*ln(N)
end

NullBIC = N*log(NullSSE/N)+log(N);
deltaBIC = NullBIC - GaussianBIC;





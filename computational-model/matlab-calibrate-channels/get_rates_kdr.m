function [alpha, beta] = get_rates_kdr(Em, params)
% Note : the starting parameters are from Wang&Buzsaki 1996

% Base model:
% alpha = 0.01 * (Em+34) / ( 1 - exp(-(Em+34)/10));

if (Em-params.alpha.E_T)==0
    alpha = params.alpha.A * params.alpha.k; % Taylor 
else
    alpha = params.alpha.A * (Em - params.alpha.E_T) / ( 1 - exp(-(Em-params.alpha.E_T)/params.alpha.k));
end

% SAME, but considered input Em as array 
% alpha = params.alpha.A * (Em - params.alpha.E_T) ./ ( 1 - exp(-(Em-params.alpha.E_T)/params.alpha.k));
% is_trap = (Em-params.alpha.E_T)==0;
% alpha(is_trap) = params.alpha.A * params.alpha.k; % Taylor 

% Base model:
% beta = 0.125 * exp(-(Em+44)/80);

beta = params.beta.B * exp(-(Em-params.beta.E_T)/params.beta.k);

end

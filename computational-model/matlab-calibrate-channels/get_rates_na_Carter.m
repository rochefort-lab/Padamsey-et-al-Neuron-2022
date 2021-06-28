function [rates_ab, rates_h] = get_Na_rates_Carter(Em, params)
% rates_ab
%   array 
%       rows: alphas (forward) / beta (backward)
%       columns: 1->K/2-1
%       3rd dim: (1) open/close (2) inactivation
% rates_h recovery/inactivation 
%   array [2 x K/2] 
%       rows: alphas (recovery) / beta (inactivation)

% K = total number of states
% K/2 = number of C/O  slash number of Inactivation states
% K/2-1 = number of transition rates alpha/beta for C/0 slash I states


subnum = 1:4; 
subnum_rev = 4:-1:1;

rates_ab = zeros(2,5,2);

alpha = params.ab.p1 *exp((Em-params.ab.Eshift)/params.ab.p2);
beta = params.ab.p3 *exp(-(Em-params.ab.Eshift)/params.ab.p4);

rates_ab(1,1:4,1) = subnum_rev*alpha;
rates_ab(2,1:4,1) = subnum*beta;

rates_ab(1,1:4,2) = subnum_rev*alpha * params.sc_a;
rates_ab(2,1:4,2) = subnum*beta / params.sc_b;

rates_ab(1,5,1:2) = params.ab.gamma;
rates_ab(2,5,1:2) = params.ab.delta;

% Recovery/Inactivation rates
rates_h = zeros(2,6);
rates_h(1,1:5) = params.h.Coff * [1, 1./(params.sc_b.^subnum)];
rates_h(2,1:5) = params.h.Con * [1, (params.sc_a.^subnum)];

rates_h(:,6) = [params.h.Ooff; params.h.Oon];

end
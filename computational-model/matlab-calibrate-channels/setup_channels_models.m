% MODEL PARAMETERS
% Standard single-compartment model for cortical neuron
%
% Model
%
% Na+NaP implementation:
%   Carter
% delayed rectifier potassium implementation:
%   H-H model for K+ from Wang & Buzsaki 1996
%
% Adjust parameters to fit data


%% UNITS

% d(-1) : c(-2) : m(-3) : u(-6) : n(-9) : p(-12)

% Conversion pS/um2 = 1e-1 mS/cm2

% Variables and units:
%   t   ms          time
%   E   mV          membrane potential
%   V   mV          displacement from resting potential
%   I   microA/cm2  current density
%   Cm  microF/cm2  membrane capacitance
%   g   mS/cm2      conductance
%   Er  mV          reversal potential


%% Setup channels - Notes about Kinetics corrections

% -- SHIFT voltage
%   E_T + v
%
% -- TAU
%   dz/dt = phi * (z_inf - z)/tau
%       => tau/phi
%

function [list_channels, list_channelnames] = setup_channels_models(gna,gkdr)


% -- Reversal potentials main ions
Er_Na = 55; % mV
Er_K = -90; % mV

%% Voltage shift for multiple channels \ Global shift to play with NaP / other channels

% Adjust the V threshold
Eshift_global = 0;
Eshift_Na = 25;
Eshift_Kdr = 15;


%% Setup Na channel - Transient Sodium

ChNa = MarkovChannel;
ChNa.fct_get_rates = @get_rates_na_Carter;
ChNa.open_state = 6;
ChNa.Er = Er_Na;
ChNa.gbar = gna; % mS/cm2 - MAX conductances

% Adjust the V threshold
Eshift = Eshift_global + Eshift_Na;
% Define the parameters to shift the kinetics as in in Wang&Buzsaki 1996
ChNa.param_rates = struct(...
    'sc_a', 2.51, 'sc_b', 5.32, ...
    'ab', struct('Eshift', Eshift, ...
    'p1', 550, 'p2', 24, ...
    'p3', 12, 'p4', 24, ...
    'gamma', 250, 'delta', 60), ...
    'h', struct(...
    'Con', 0.01, 'Coff', 2, ...
    'Oon', 8, 'Ooff', 0.05) );


%% Setup Kdr channel - Delayed Rectifier Potassium


ChKdr = BasicChannel;
ChKdr.fct_get_rates = @get_rates_kdr;
ChKdr.fct_get_open_proba = @(z)z^4;
ChKdr.Er = Er_K;
ChKdr.gbar = gkdr; % mS/cm2 - MAX conductances 

% Adjust the V threshold
Eshift = Eshift_global + Eshift_Kdr;

% Define the parameters - adjust the kinetics as in in Wang&Buzsaki 1996
Phi = 5;
ChKdr.param_rates = struct(...
    'alpha', struct(...
    'A', 0.01*Phi, 'E_T', -34+Eshift, 'k', 10), ...
    'beta', struct(...
    'B', 0.125*Phi, 'E_T', -44+Eshift, 'k', 80));


%% Save the 2 channels

list_channels = {ChNa, ChKdr};
list_channelnames = {'Na', 'Kdr'};


end


% Variables and units:
%   t   ms          time
%   Em  mV          membrane potential
%   I   microA/cm2  current density
%   Cm  microF/cm2  membrane capacitance
%   g   mS/cm2      conductance
%   Er  mV          reversal potential


function [Y, t, Z_track, IZ_track] = run_trace(dt, T_simu, pas_params, list_channels, Em0, input_inject)


% -- List of channels
% Total number of currents included in membrane current
% Nchannels = length(list_channels);
% Reshape to make sure we can use cellfun combined with Z
if size(list_channels,1)>1
    list_channels = list_channels(:)';
end

% -- Initialise convenient variables

% Concatenate the reversing / equilibrium potential
Eeq = cellfun(@(ch)ch.Er, list_channels, 'un', true); 
% [1 x Ncurrents]

Gbar = cellfun(@(ch)ch.gbar, list_channels, 'un', true);
% [1 x Ncurrents]

% -- Initialise the output variables --

% Number of time points
Nt = length(input_inject)+1;
t = 0:dt:T_simu;

Y = zeros(Nt, 3); % membrane potential, current, conductance
% Artificially add a datapoint to store the initial values. Corresponds to
% t=0

% Initialise membrane potential to given Em0
Y(1,1) = Em0;

% Compute initial conditions
[Z,~,~,g0,ImZ0] = cellfun(@(ch)ch.compute_all_inf(Em0), list_channels, 'un', false);

% Membrane currents 
Y(1,2) = pas_params.gL*( Em0 - pas_params.ErL ) + sum([ImZ0{:}]);

IZ_track = zeros(Nt,2);
IZ_track(1,:) = [ImZ0{:}];

% Membrane conductance
GZ = [g0{:}];
Y(1,3) = pas_params.gL + sum(GZ);

n_states_z = cellfun(@(z)size(z,2), Z, 'un', true);
n_states_tot = sum(n_states_z);
Z_track = zeros(Nt,n_states_tot);
Z_track(1,:)=[Z{:}];


% -- Run --

for k_t=2:Nt
    

    % Update membrane potential
    E_inf = (GZ*Eeq' + pas_params.gL*pas_params.ErL + input_inject(k_t-1)) / (sum(GZ)+pas_params.gL);
    tau_v = pas_params.Cm / (sum(GZ)+pas_params.gL);
    Em = E_inf + (Y(k_t-1,1) - E_inf)*exp(-dt/tau_v);
    
    Y(k_t,1) = Em;
    
    % TO DO: membrane current!
    
    % Update gating variables t, Em, z
    Z = cellfun(@(ch,x)ch.update_state(dt, Em, x), list_channels, Z, 'un', false);
    % Update the open state
    Popen = cellfun(@(ch,x)ch.get_open_proba(x), list_channels, Z, 'un', true);
    % Update the conductances
    GZ = Popen.*Gbar;

    Y(k_t,3) = pas_params.gL + sum(GZ);
    
    ImZ = cellfun(@(ch,ginf)ginf*(Em-ch.Er), list_channels, num2cell(GZ), 'un', false);
        
    Z_track(k_t,:)=[Z{:}];
    IZ_track(k_t,:) = [ImZ{:}];
    
end

Z_track = mat2cell(Z_track, Nt, n_states_z);


end




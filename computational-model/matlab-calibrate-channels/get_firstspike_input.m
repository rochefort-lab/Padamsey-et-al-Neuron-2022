% Variables and units:
%   t   ms          time
%   Em  mV          membrane potential
%   I   microA/cm2  current density
%   Cm  microF/cm2  membrane capacitance
%   g   mS/cm2      conductance
%   Er  mV          reversal potential


function [onInput, offInput, is_success] = get_firstspike_input(dt, pas_params, list_channels, Em0, input_inject0, APthreshold, precision, onInput, offInput)

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

% Number of time points
Nt = length(input_inject0);

current_input = round(((offInput+onInput)/2)/(precision/10))*precision/10;
while (onInput-offInput)>precision
    
    input_inject = current_input * input_inject0;

    % Init
    Em=Em0;
    [Z,~,~,g0,~] = cellfun(@(ch)ch.compute_all_inf(Em0), list_channels, 'un', false);
    % Init the conductance
    GZ = [g0{:}];

    k_t = 0;
    has_ap = false;
    % Monitor membrane potential until reach desired threshold
    while k_t<Nt && ~has_ap
        
        k_t = k_t+1;
        [Em, Z, GZ] = inject_current(dt, pas_params, list_channels, Eeq, Gbar, Z, GZ, Em, input_inject(k_t));
        has_ap = Em>=APthreshold;
        
    end
    % Check if has generated an AP
    if ~has_ap
        % If not AP, we need to increase the input
        offInput = current_input;
        current_input = round(((current_input+onInput)/2)/(precision/10))*precision/10;
    else
        % If AP, we need to decrease the input
        onInput = current_input;
        current_input = round(((current_input+offInput)/2)/(precision/10))*precision/10;
    end
end



% Check if the on is on...
input_inject = onInput * input_inject0;    
Em=Em0;
[Z,~,~,g0,~] = cellfun(@(ch)ch.compute_all_inf(Em0), list_channels, 'un', false);
% Init the conductance
GZ = [g0{:}];
k_t = 0;
has_ap = false;
while k_t<Nt && ~has_ap
    k_t = k_t+1;
    [Em, Z, GZ] = inject_current(dt, pas_params, list_channels, Eeq, Gbar, Z, GZ, Em, input_inject(k_t));
    has_ap = Em>=APthreshold; 
end
n1 = has_ap;
%  ... and off is off 
input_inject = offInput * input_inject0;
Em=Em0;
[Z,~,~,g0,~] = cellfun(@(ch)ch.compute_all_inf(Em0), list_channels, 'un', false);
% Init the conductance
GZ = [g0{:}];
k_t = 0;
has_ap = false;
while k_t<Nt && ~has_ap
    k_t = k_t+1;
    [Em, Z, GZ] = inject_current(dt, pas_params, list_channels, Eeq, Gbar, Z, GZ, Em, input_inject(k_t));
    has_ap = Em>=APthreshold; 
end
n0 = ~has_ap;

is_success = n1 && n0;




end

function [Em, Z, GZ] = inject_current(dt, pas_params, list_channels, Eeq, Gbar, Z, GZ, Em0, Iinject)

% Update membrane potential
E_inf = (GZ*Eeq' + pas_params.gL*pas_params.ErL + Iinject) / (sum(GZ)+pas_params.gL);
tau_v = pas_params.Cm / (sum(GZ)+pas_params.gL);
Em = E_inf + (Em0 - E_inf)*exp(-dt/tau_v);


% Update gating variables t, Em, z
Z = cellfun(@(ch,x)ch.update_state(dt, Em, x), list_channels, Z, 'un', false);
% Update the open state
Popen = cellfun(@(ch,x)ch.get_open_proba(x), list_channels, Z, 'un', true);
% Update the conductances
GZ = Popen.*Gbar;



end



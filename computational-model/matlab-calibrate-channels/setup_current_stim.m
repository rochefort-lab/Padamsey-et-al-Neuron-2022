function [input_inject0] = setup_current_stim(dt, T_simu, stim_params)

% Number of time points
Nt = round(T_simu/dt);

% -- Initialise the stimulus: inject current step or custom waveform --

% Setup the injected current
if strcmp(stim_params.type, 'i_step')
    % Current step
    input_inject0 = zeros(1,Nt);
    if ~isempty(stim_params.t_Ie_inj) 
        k0 = max(round(stim_params.t_Ie_inj(1)/dt), 1);
        k1 = min(round(stim_params.t_Ie_inj(2)/dt), Nt);
        input_inject0(k0:k1) = 1;
    end
else
    t_i = dt:dt:T_simu;
    input_inject0 = setupdoubleexp(t_i, stim_params.tau_decay, stim_params.tau_rise);
end



end


function [x] = setupdoubleexp(t_i, tau_decay, tau_rise)

x = exp(-t_i/tau_decay) - exp(-t_i/tau_rise);

x = x/max(x);

end
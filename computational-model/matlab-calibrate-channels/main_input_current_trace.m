format long


%% USER INPUT: Setup parameters

dt = 0.01;
T_simu=150;

% -- Model

model_params = struct();

model_params.name = 'control';
model_params.color = [0.4,0.4,0.4];
model_params.pas_params = struct('gL', 0.1, 'ErL', -75, 'Cm', 1); % basic: -65
model_params.Em0 = -75; 

% --- Channels 
gna = 26;
gkdr = 6 ;
[list_channels, list_channelnames] = setup_channels_models(gna,gkdr);


% --- Stim input current:

% * double exp
% stim_params = struct('type', 'isyn');
% stim_params.tau_rise = 50 ; %50;
% stim_params.tau_decay = 100 ; %100; % in ms

% * step
stim_params = struct('type', 'i_step');
stim_params.t_Ie_inj = [5,120]; % in ms


%% Setup stim

[input_inject0] = setup_current_stim(dt, T_simu, stim_params);


%% Run Trace

s = 2.4;


[Y, t, Z_track, IZ_track] = run_trace(dt, T_simu, ...
    model_params.pas_params, list_channels, ...
    model_params.Em0, s*input_inject0);


%%

lw = 1.5;

figure(100)
clf
subplot(3,2,1)
plot(t,[0,s*input_inject0], 'LineWidth', lw, 'Color', model_params.color)
title('Input current')

subplot(3,2,3)
plot(t,Y(:,1), 'LineWidth', lw, 'Color', model_params.color)
% ylim([model_params.control.Em0-5, 20])
ylim([model_params.Em0-5, 50])
title('Membrane potential (mV)')

subplot(3,2,4)
plot(t,Y(:,3), 'LineWidth', lw, 'Color', model_params.color)
title('Membrane conductance')

subplot(3,2,6)
hold on
plot(t,Z_track{2}, 'LineWidth', lw, 'Color', 'b')
plot(t,Z_track{1}(:,list_channels{1}.open_state), 'LineWidth', lw, 'Color', 'r')
legend({'Kdr (z)', 'Na (Open state)'})
title('Channels states')

subplot(3,2,5)
hold on
plot(t,-IZ_track(:,1), 'LineWidth', lw, 'Color', 'r')
plot(t,IZ_track(:,2), 'LineWidth', lw, 'Color', 'b')
legend({'Kdr', '-Na'})
title('Channels currents')


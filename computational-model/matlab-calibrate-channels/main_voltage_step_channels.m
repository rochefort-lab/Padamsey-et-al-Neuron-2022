format long


%% USER INPUT: Setup parameters

E_steps = -90:1:20;


% --- Channels 
gna = 26;
gkdr = 6 ;
[list_channels, list_channelnames] = setup_channels_models(gna,gkdr);



%% Steps to get activation curve

n = length(E_steps);

Z = zeros(n,2);
Tau = zeros(n,2);

for k=1:n
    
    u = cellfun(@(ch,x)ch.compute_steadystate(E_steps(k)), list_channels, 'un', false);
    
    Z(k,1) = u{1}(:,list_channels{1}.open_state);
    Z(k,2) = u{2};
    
    Tau(k,2) = list_channels{2}.compute_tau(E_steps(k));
    
end

Z(:,1)= Z(:,1)/max(Z(:,1));

%%

lw = 2;

figure(400)
clf
subplot(2,1,1)
plot(E_steps,Z, 'LineWidth', lw)
title('Steady state')
legend({'Na (Persistent Po, norm.)', 'Kdr'}, 'location', 'northwest')

subplot(2,1,2)
plot(E_steps,Tau, 'LineWidth', lw)
title('Tau (Kdr only)')
legend({'Na (none)', 'Kdr'})



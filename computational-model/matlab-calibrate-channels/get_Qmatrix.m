function [Q] = get_Qmatrix(rates_ab, rates_h) % precision

% ATTENTION for singular matrix, test using rank
% see: https://stackoverflow.com/questions/13145948/how-to-find-out-if-a-matrix-is-singular

% ATTENTION - double numerical incaccuracies!!

% rates_ab
%   array 
%       rows: alphas (forward) / beta (backward)
%       columns: 1->K/2-1
%       3rd dim: (1) open/close (2) inactivation
% rates_h recovery/inac 
%   array [2 x K/2] 
%       rows: alphas (forward) / beta (backward)

% K = total number of states
% K/2 = number of C/O  slash number of Inactivation states
% K/2-1 = number of transition rates alpha/beta for C/0 slash I states


Ntr = size(rates_ab,2); % number of transitions
Nshalf = Ntr+1;% half number of states

% Close-Close
% Directly using matlab diag function...
CO1 = diag(rates_ab(1,:,1),1);
CO2 = diag(rates_ab(2,:,1),-1);
% CO1 = [[zeros(Ntr,1),diag(rates_ab(1,:,1))] ; zeros(1,Nshalf)];
% CO2 = [zeros(1,Nshalf);[diag(rates_ab(2,:,1)),zeros(Ntr,1)]];

if isempty(rates_h)
    Q = CO1+CO2;
else
    % Directly using matlab diag function...
    I1 = diag(rates_ab(1,:,2),1);
    I2 = diag(rates_ab(2,:,2),-1);
%     I1 = [[zeros(Ntr,1),diag(rates_ab(1,:,2))] ; zeros(1,Nshalf)];
%     I2 = [zeros(1,Nshalf);[diag(rates_ab(2,:,2)),zeros(Ntr,1)]];

    Q = [ [CO1+CO2 , diag(rates_h(2,:))] ; [diag(rates_h(1,:)) , I1+I2]];
end


% Now set the elements in the diagonal so that sum each row =0
diag_q_el = -sum(Q,2);

Q = Q + diag(diag_q_el);

% % See what is faster : above or below
% Q(logical(eye(Khalf*2)))=-sum(Q,2);

end

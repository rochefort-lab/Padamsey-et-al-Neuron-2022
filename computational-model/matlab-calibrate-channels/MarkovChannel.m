classdef MarkovChannel
    
    properties
        
        fct_get_rates
        param_rates
        open_state {mustBeNumeric}
        gbar
        Er
        
    end
    
    methods
        
        
        function [popen] = get_open_proba(obj, z)
            
            popen = z(obj.open_state);
            
        end
        

        function [zinf, tauz, pinf, ginf, iinf] = compute_all_inf(obj, Em)
            
            % State of gating vars 
            zinf = obj.compute_steadystate(Em);
            
            % Add dummy tau
            tauz = zeros(size(zinf));
            
            % Open probability
            pinf = obj.get_open_proba(zinf);
            % Conductance
            ginf = obj.gbar * pinf;
            % Membrane current - voltage dependent conductances
            iinf = ginf * (Em - obj.Er);
            
        end
        
        function [zinf] = compute_steadystate(obj, Em)
            
            % Get rates
            [rates_ab, rates_h] = obj.fct_get_rates(Em, obj.param_rates);
            
            % Get Q-matrix
            % See Colqhoun&Hawkes
            Q = get_Qmatrix(rates_ab, rates_h);

            K = size(Q,1);
            S = [Q, ones(K,1)];
            
            zinf = ones(1,K) / (S*S');
            
        end
        
        function [z] = update_state(obj, t, Em, z)
            % p0 : [1 x K]
            
            % Get rates
            [rates_ab, rates_h] = obj.fct_get_rates(Em, obj.param_rates);
            
            % Get Q-matrix
            % See Colqhoun&Hawkes
            Q = get_Qmatrix(rates_ab, rates_h);

            % Get the eigenvectors / eigenvalues
            [X,lambda] = eig(-Q);
            % [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V
            % whose columns are the corresponding right eigenvectors, so that A*V =
            % V*D.
            lambda = exp(-t*diag(lambda));
            
            % Get the spectral matrices
            Y = X^(-1);
            K=length(z);
            A = zeros(K,K,K);
            for r=1:K
                A(:,:,r) = lambda(r) * X(:,r)*Y(r,:);
            end
            
            z = z*sum(A,3);
            
        end
        
        
    end
    
end
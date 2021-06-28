classdef BasicChannel
    
    properties
        
        fct_get_rates
        fct_get_open_proba
        param_rates
        gbar
        Er
        
    end
    
    methods
        
        
        function [popen] = get_open_proba(obj, z)
            
            popen = obj.fct_get_open_proba(z);
            
        end
        
        function [zinf, tauz, pinf, ginf, iinf] = compute_all_inf(obj, Em)
            % State of gating vars
            [zinf, tauz] = obj.compute_steadystate_tau(Em);
            % Open probability
            pinf = obj.get_open_proba(zinf);
            % Conductance
            ginf = obj.gbar * pinf;
            % Membrane current - voltage dependent conductances
            iinf = ginf * (Em - obj.Er);
            
        end
        
        function [zinf, tauz] = compute_steadystate_tau(obj, Em)
            
            [alpha, beta] = obj.fct_get_rates(Em, obj.param_rates);
            
            [zinf] = obj.compute_steadystate(alpha, beta);
            
            [tauz] = obj.compute_tau(alpha, beta);
            
        end
        
        function [zinf] = compute_steadystate(obj, varargin)
            if nargin==2
                [alpha, beta] = obj.fct_get_rates(varargin{1}, obj.param_rates);
            else
                alpha=varargin{1}; beta=varargin{2};
            end
            
            zinf = alpha ./ (alpha + beta);
        
        end
        
        function [tauz] = compute_tau(obj, varargin)
            if nargin==2
                [alpha, beta] = obj.fct_get_rates(varargin{1}, obj.param_rates);
            else
                alpha=varargin{1}; beta=varargin{2};
            end
            
            tauz = 1 ./ (alpha + beta);

        end
        
        function [z] = update_state(obj, t, Em, z)
            
            % Get steady state and tau
            [zinf, tauz] = obj.compute_steadystate_tau(Em);
            
            z = zinf + (z-zinf).*exp(-t./tauz) ; % [1 x nvariables]
            
            % Backup - case with tau = 0 instantaneous
            
%             is_instant = isnan(tauz);
%             if any(is_instant)
%                 z(~is_instant) = zinf(~is_instant) + (z(~is_instant)-zinf(~is_instant)).*exp(-t./tauz(~is_instant)) ; % [1 x nvariables]
%                 z(is_instant) = zinf(is_instant) ; % [1 x nvariables]
%             else
%                 z = zinf + (z-zinf).*exp(-t./tauz) ; % [1 x nvariables]
%             end
            
        end
        
        
    end
    
end
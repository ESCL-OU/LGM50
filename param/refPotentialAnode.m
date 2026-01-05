%% Reference Potential for Pos. Electrode: Unref(theta_n)
%   Created May 12, 2025 by Wonoo Choo

function [Uref,varargout] = refPotentialAnode(p,theta)

% % Polynomail Fit
% Uref = ppvalFast(p.Uppp,theta);

% DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
Uref = 1.9793 * exp(-39.3631 * theta) ...
     + 0.2482 ...
     - 0.09069 * tanh(29.8538 * (theta - 0.1234))  ...
     - 0.04478 * tanh(14.9159 * (theta - 0.2769)) ...
     - 0.0205 * tanh(30.4444 * (theta - 0.6103));

% Gradient of OCP wrt theta
if(nargout == 2)
    
%     % Polynomail Fit
%     dUref = ppvalFast(p.dUppp,theta);
%     varargout{1} = dUref / p.c_s_p_max;

dUref = -77.9176 * exp(-39.3631 * theta) ...
        - 2.70744 * sech(29.8538 * (theta - 0.1234)).^2 ...
        - 0.6679 * sech(14.9159 * (theta - 0.2769)).^2 ...
        - 0.6241 * sech(30.4444 * (theta - 0.6103)).^2;

varargout{1} = dUref;
    
end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end
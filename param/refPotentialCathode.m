%% Reference Potential for Pos. Electrode: Unref(theta_n)
%   Created May 12, 2025 by Wonoo Choo

function [Uref,varargout] = refPotentialCathode(p,theta)

Uref = -0.8090 * theta ...
       + 4.4875 ...
       - 0.0428 * tanh(18.5138 * (theta - 0.5542)) ...
       - 17.7326 * tanh(15.7890 * (theta - 0.3117)) ...
       + 17.5842 * tanh(15.9308 * (theta - 0.3120));
% Gradient of OCP wrt theta
if(nargout == 2)
    
dUref = -0.8090 ...
        - 0.7924 * sech(18.5138 * (theta - 0.5542)).^2 ...
        - 279.98  * sech(15.7890 * (theta - 0.3117)).^2 ...
        + 280.13  * sech(15.9308 * (theta - 0.3120)).^2;
varargout{1} = dUref;
    
end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end
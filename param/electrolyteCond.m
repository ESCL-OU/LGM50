%% Electrolyte Conductivity Function: kappa(c_e) [1/Ohms*m]
%   Created May 15, 2025 by Wonoo Choo

function [kappa,varargout] = electrolyteCond(p, c_e, T)
%  Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
%     comes from [1].

%     References
%     ----------
%     .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
%     modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
%     Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.

arrh_k = exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T));

kappa = 0.1297 * (c_e / 1000) .^ 3 - 2.51 * (c_e / 1000) .^ 1.5 + 3.329 * (c_e / 1000);
kappa = kappa * arrh_k;

% Nyman et al. (2008) does not provide temperature dependence

if(nargout == 2)
    dkappa =  0.1297 * 3 * (c_e / 1000)^2 / 1000 ...
       - 2.51 * 1.5 * (c_e / 1000).^0.5 / 1000 ...
       + 3.329 / 1000;
    varargout{1} = dkappa * arrh_k;
end
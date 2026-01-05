%% Electrolyte Diffusion Coefficient Function: D_e(c_e) [m^2/s]
%   Created May 15, 2025 by Wonoo Choo

function [D_e,varargout] = electrolyteDe(p, c_e, T)
% Diffusivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
% comes from [1], with Arrhenius temperature dependence added from [2].

% References
% ----------
% .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
% modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
% Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.
% .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
% a lithium-ion battery i. determination of parameters." Journal of the
% Electrochemical Society 162.9 (2015): A1836-A1848.

arrh_De = exp(p.E.De/p.R*(1/p.T_ref - 1/T));

D_e = 8.794e-11 * (c_e / 1000) .^ 2 - 3.972e-10 * (c_e / 1000) + 4.862e-10;
D_e = D_e * arrh_De;

if(nargout == 2)
    dD_e = 2 * 8.794e-11 * (c_e / 1000) / 1000 - 3.972e-10 / 1000;
    varargout{1} = dD_e * arrh_De;
end
%% Exchange Current Density function i_0
%   Created May 14, 2025 by Wonoo Choo

function [i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e,T)

% Parse out concentrations in anode and cathode
c_e_n = c_e(1:p.Nxn-1);
c_e_p = c_e(p.Nxn+p.Nxs-1:end);

% Compute exchange current density
i_0n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T)) * ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^p.alph;
i_0p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T)) * ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_e_p).^p.alph;
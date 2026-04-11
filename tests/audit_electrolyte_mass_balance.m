function diagout = audit_electrolyte_mass_balance(c_e, jn, jp, p)
%AUDIT_ELECTROLYTE_MASS_BALANCE
% Diagnose where discrete electrolyte-lithium non-conservation comes from
% in the current semi-discrete ODE used by ode_electrolyte.m.
%
% This does NOT merely check whether conservation is violated.
% It decomposes dNe_disc/dt = w' * c_e_dot into contributions from:
%   1) dD_eff/dc * (dc/dx)^2
%   2) D_eff * d2c/dx2
%   3) source terms from jn / jp
%
% Inputs
%   c_e : electrolyte interior state vector, size (p.Nx-3)x1
%   jn  : anode interfacial molar flux used in ode_electrolyte
%   jp  : cathode interfacial molar flux used in ode_electrolyte
%   p   : parameter struct
%
% Output fields
%   .dNe_total
%   .regions.anode / .sep / .cathode
%   .terms.*  (global sums)
%   .vectors.* (nodewise contribution vectors)
%
% Notes
% - This matches the implementation in ode_electrolyte.m exactly.
% - It uses only the interior state vector c_e, because that is what the ODE uses.

    c_e = c_e(:);

    [M1n,M2n,M3n,M4n,M5n, ...
     M1s,M2s,M3s,M4s, ...
     M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

    % Boundary/interface concentrations reconstructed algebraically
    c_e_bcs = C_ce * c_e;

    % Split by region, same as ode_electrolyte.m
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;

    idx_n = 1:Nn;
    idx_s = (Nn+1):(Nn+Ns);
    idx_p = (Nn+Ns+1):(Nn+Ns+Np);

    c_en = c_e(idx_n);
    c_es = c_e(idx_s);
    c_ep = c_e(idx_p);

    % Effective diffusivity and derivative, same as ode_electrolyte.m
    [D_en,dD_en] = electrolyteDe(p, c_en, p.T_amb);
    [D_es,dD_es] = electrolyteDe(p, c_es, p.T_amb);
    [D_ep,dD_ep] = electrolyteDe(p, c_ep, p.T_amb);

    D_en_eff  = D_en  .* p.epsilon_e_n.^(p.brug-1);
    dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

    D_es_eff  = D_es  .* p.epsilon_e_s.^(p.brug-1);
    dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

    D_ep_eff  = D_ep  .* p.epsilon_e_p.^(p.brug-1);
    dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

    % --- Reconstruct the exact three pieces in ode_electrolyte.m ---

    % First-derivative approximations
    gx_n = M1n*c_en + M2n*c_e_bcs(1:2);
    gx_s = M1s*c_es + M2s*c_e_bcs(2:3);
    gx_p = M1p*c_ep + M2p*c_e_bcs(3:4);

    % Second-derivative approximations
    gxx_n = M3n*c_en + M4n*c_e_bcs(1:2);
    gxx_s = M3s*c_es + M4s*c_e_bcs(2:3);
    gxx_p = M3p*c_ep + M4p*c_e_bcs(3:4);

    % Term 1: dD_eff/dc * (dc/dx)^2
    term1_n = dD_en_eff .* (gx_n.^2);
    term1_s = dD_es_eff .* (gx_s.^2);
    term1_p = dD_ep_eff .* (gx_p.^2);

    % Term 2: D_eff * d2c/dx2
    term2_n = D_en_eff .* gxx_n;
    term2_s = D_es_eff .* gxx_s;
    term2_p = D_ep_eff .* gxx_p;

    % Term 3: source terms
    src_n = diag(M5n) * jn;
    src_s = zeros(Ns,1);
    src_p = diag(M5p) * jp;

    % Full c_e_dot reconstructed exactly like ode_electrolyte.m
    c_en_dot = term1_n + term2_n + src_n;
    c_es_dot = term1_s + term2_s + src_s;
    c_ep_dot = term1_p + term2_p + src_p;
    c_e_dot  = [c_en_dot; c_es_dot; c_ep_dot];

    % Weighted discrete lithium content based on INTERIOR states
    Dx_n = p.L_n * p.delta_x_n;
    Dx_s = p.L_s * p.delta_x_s;
    Dx_p = p.L_p * p.delta_x_p;

    w_n = p.epsilon_e_n * Dx_n * ones(Nn,1);
    w_s = p.epsilon_e_s * Dx_s * ones(Ns,1);
    w_p = p.epsilon_e_p * Dx_p * ones(Np,1);
    w   = [w_n; w_s; w_p];

    % Global derivative
    dNe_total = w.' * c_e_dot;

    % Decompose by term
    vec_term1 = [term1_n; term1_s; term1_p];
    vec_term2 = [term2_n; term2_s; term2_p];
    vec_src   = [src_n;   src_s;   src_p];

    dNe_term1 = w.' * vec_term1;
    dNe_term2 = w.' * vec_term2;
    dNe_src   = w.' * vec_src;

    % Decompose by region
    dNe_n = w_n.' * c_en_dot;
    dNe_s = w_s.' * c_es_dot;
    dNe_p = w_p.' * c_ep_dot;

    % More detailed region-by-term decomposition
    reg.anode.term1   = w_n.' * term1_n;
    reg.anode.term2   = w_n.' * term2_n;
    reg.anode.source  = w_n.' * src_n;
    reg.anode.total   = dNe_n;

    reg.sep.term1     = w_s.' * term1_s;
    reg.sep.term2     = w_s.' * term2_s;
    reg.sep.source    = w_s.' * src_s;
    reg.sep.total     = dNe_s;

    reg.cathode.term1 = w_p.' * term1_p;
    reg.cathode.term2 = w_p.' * term2_p;
    reg.cathode.source= w_p.' * src_p;
    reg.cathode.total = dNe_p;

    % Package outputs
    diagout = struct();
    diagout.dNe_total = dNe_total;

    diagout.terms = struct( ...
        'term1_dD_grad2', dNe_term1, ...
        'term2_D_lap',    dNe_term2, ...
        'source',         dNe_src);

    diagout.regions = reg;

    diagout.vectors = struct( ...
        'w', w, ...
        'c_e_dot', c_e_dot, ...
        'term1', vec_term1, ...
        'term2', vec_term2, ...
        'source', vec_src, ...
        'gx_n', gx_n, 'gx_s', gx_s, 'gx_p', gx_p, ...
        'gxx_n', gxx_n, 'gxx_s', gxx_s, 'gxx_p', gxx_p);

    diagout.boundary_values = c_e_bcs;
end
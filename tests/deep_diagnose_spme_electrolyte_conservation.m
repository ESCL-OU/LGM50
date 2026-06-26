function report = deep_diagnose_spme_electrolyte_conservation()
%DEEP_DIAGNOSE_SPME_ELECTROLYTE_CONSERVATION
% Diagnose whether electrolyte RHS conserves electrolyte-phase Li inventory.
%
% It tests:
%   1) diffusion-only drift at I = 0
%   2) source-only drift from a uniform electrolyte profile
%   3) full RHS drift from a perturbed electrolyte profile
%
% Run from repo root or from tests folder after adding paths.
%
% Interpretation:
%   dNdt_state = sum_j eps_j * dx_j * sum(cdot_j)
%   For a conservative semi-discrete electrolyte model, dNdt_state should be
%   near zero for no side reactions, up to solver/roundoff.

    base_dir = fileparts(mfilename('fullpath'));
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    init_soc = 0.5;
    t_amb = 25;
    Nr = 15;
    delta_t = 1;
    Nxn = 20;
    Nxs = 5;
    Nxp = 20;
    disable_thermal = true;

    p = load_SPMe_params(init_soc, t_amb, Nr, delta_t, Nxn, Nxs, Nxp, ...
        disable_thermal, 'none', 'fdm', false, false, false);

    I = 0.2 * p.OneC; % A, same sign convention as SPMe.m
    jn = I / (p.Faraday * p.a_s_n * p.Area_n * p.L_n);
    jp = -I / (p.Faraday * p.a_s_p * p.Area_p * p.L_p);

    Nx = p.Nx - 3;
    ce_uniform = p.c_e * ones(Nx, 1);

    rng(1);
    ce_perturbed = ce_uniform .* (1 + 0.02 * randn(Nx,1));
    ce_perturbed = max(1, ce_perturbed);

    report = struct();

    report.fdm.diffusion_only = eval_case(@ode_electrolyte, ce_perturbed, 0, 0, p);
    report.fvm.diffusion_only = eval_case(@ode_electrolyte_conservative, ce_perturbed, 0, 0, p);

    report.fdm.source_only = eval_case(@ode_electrolyte, ce_uniform, jn, jp, p);
    report.fvm.source_only = eval_case(@ode_electrolyte_conservative, ce_uniform, jn, jp, p);

    report.fdm.full_rhs = eval_case(@ode_electrolyte, ce_perturbed, jn, jp, p);
    report.fvm.full_rhs = eval_case(@ode_electrolyte_conservative, ce_perturbed, jn, jp, p);

    Nn = p.Nxn - 1; Np = p.Nxp - 1;
    dxn = p.L_n * p.delta_x_n;
    dxp = p.L_p * p.delta_x_p;

    report.source_balance_state = ...
        Nn * dxn * (1 - p.t_plus) * p.a_s_n * jn + ...
        Np * dxp * (1 - p.t_plus) * p.a_s_p * jp;

    report.source_balance_physical = ...
        p.L_n * (1 - p.t_plus) * p.a_s_n * jn + ...
        p.L_p * (1 - p.t_plus) * p.a_s_p * jp;

    disp(report)
end

function s = eval_case(rhsfun, ce, jn, jp, p)
    cdot = rhsfun(ce, jn, jp, p);
    cdot = full(double(cdot));

    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;

    cdot_n = cdot(1:Nn);
    cdot_s = cdot(Nn+1:Nn+Ns);
    cdot_p = cdot(Nn+Ns+1:Nn+Ns+Np);

    dxn = p.L_n * p.delta_x_n;
    dxs = p.L_s * p.delta_x_s;
    dxp = p.L_p * p.delta_x_p;

    s = struct();
    s.dNdt_state = ...
        p.epsilon_e_n * dxn * sum(cdot_n) + ...
        p.epsilon_e_s * dxs * sum(cdot_s) + ...
        p.epsilon_e_p * dxp * sum(cdot_p);
    s.region_terms = [ ...
        p.epsilon_e_n * dxn * sum(cdot_n), ...
        p.epsilon_e_s * dxs * sum(cdot_s), ...
        p.epsilon_e_p * dxp * sum(cdot_p)];
    s.max_abs_cdot = max(abs(cdot));
end

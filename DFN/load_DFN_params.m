function p = load_DFN_params(init_SOC, t_amb, delta_t, disable_thermal, varargin)

%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LGM50

if ~isscalar(init_SOC) || init_SOC < 0 || init_SOC > 1
    error('init_SOC must be a scalar in [0, 1].');
end

% Optional discretization/degradation override:
% load_DFN_params(..., Nr, Nxn, Nxs, Nxp, degradation_model,
%                 enable_lithium_plating, enable_anode_porosity_loss)
Nr_opt = p.Nr;
Nxn_opt = p.Nxn;
Nxs_opt = p.Nxs;
Nxp_opt = p.Nxp;
if numel(varargin) >= 4
    if ~isempty(varargin{1}), Nr_opt = varargin{1}; end
    if ~isempty(varargin{2}), Nxn_opt = varargin{2}; end
    if ~isempty(varargin{3}), Nxs_opt = varargin{3}; end
    if ~isempty(varargin{4}), Nxp_opt = varargin{4}; end
end
if numel(varargin) >= 5 && ~isempty(varargin{5})
    p.degradation_model = lower(char(varargin{5}));
else
    p.degradation_model = 'none';
end
if ~ismember(p.degradation_model, {'none', 'safari'})
    error('DFN degradation_model must be one of: none, safari.');
end
if numel(varargin) >= 6 && ~isempty(varargin{6})
    p.enable_lithium_plating = logical(varargin{6});
else
    p.enable_lithium_plating = false;
end
if numel(varargin) >= 7 && ~isempty(varargin{7})
    p.enable_anode_porosity_loss = logical(varargin{7});
else
    p.enable_anode_porosity_loss = false;
end
if numel(varargin) >= 8 && ~isempty(varargin{8})
    p.dfn_lpl_coupling_mode = lower(char(varargin{8}));
else
    p.dfn_lpl_coupling_mode = 'spme_compatible';
end
if ~ismember(p.dfn_lpl_coupling_mode, {'spme_compatible', 'yang_total_flux'})
    error('p.dfn_lpl_coupling_mode must be one of: spme_compatible, yang_total_flux.');
end
p.dfn_degradation_enabled = ~strcmpi(char(p.degradation_model), 'none');
p.enable_lithium_plating = p.dfn_degradation_enabled && p.enable_lithium_plating;
p.enable_anode_porosity_loss = p.dfn_degradation_enabled && p.enable_anode_porosity_loss;
if ~isfield(p, 'kappa_SEI') || isempty(p.kappa_SEI)
    p.kappa_SEI = p.kappa_sei;
end
if ~isfield(p, 'delta_sei_min') || isempty(p.delta_sei_min)
    p.delta_sei_min = 1e-12;
end
if ~isfield(p, 'epsilon_sei') || isempty(p.epsilon_sei)
    p.epsilon_sei = 0.05;
end
if ~isfield(p, 'epsilon_e_n_min') || isempty(p.epsilon_e_n_min)
    p.epsilon_e_n_min = 0.03;
end
if ~isfield(p, 'c_e_min') || isempty(p.c_e_min)
    p.c_e_min = 1;
end
if ~isfield(p, 'dfn_exp_arg_min') || isempty(p.dfn_exp_arg_min)
    p.dfn_exp_arg_min = -50;
end
if ~isfield(p, 'dfn_exp_arg_max') || isempty(p.dfn_exp_arg_max)
    p.dfn_exp_arg_max = 50;
end
if ~isfield(p, 'dfn_lpl_warn_flux') || isempty(p.dfn_lpl_warn_flux)
    p.dfn_lpl_warn_flux = 1e-4;
end
if ~isfield(p, 'dfn_lpl_max_iter') || isempty(p.dfn_lpl_max_iter)
    p.dfn_lpl_max_iter = 25;
end
if ~isfield(p, 'dfn_lpl_tol') || isempty(p.dfn_lpl_tol)
    p.dfn_lpl_tol = 1e-12;
end

% DFN degradation Jacobian mode.
% 'finite_difference' is safest but very slow.
% 'padded_analytic' uses Scott's original analytic DFN Jacobian and freezes
% degradation derivatives in the Newton Jacobian. The full DAE residual still
% contains SEI, plating, porosity, and electrolyte feedback.
if numel(varargin) >= 9 && ~isempty(varargin{9})
    p.dfn_jacobian_mode = lower(char(varargin{9}));
else
    p.dfn_jacobian_mode = 'padded_analytic';
end

if ~ismember(p.dfn_jacobian_mode, {'finite_difference', 'padded_analytic'})
    error('p.dfn_jacobian_mode must be one of: finite_difference, padded_analytic.');
end

%%% INITIAL CONDITIONS
Nr = Nr_opt;
Nxn = Nxn_opt;
Nxs = Nxs_opt;
Nxp = Nxp_opt;

% Apply discretization first so all downstream matrices/states are consistent.
p.Nr = Nr;
p.Nxn = Nxn;
p.Nxs = Nxs;
p.Nxp = Nxp;
p.Nx = p.Nxn+p.Nxs+p.Nxp;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

% Initialize from intended SoC via SoC->Veq->init_cs.
repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_dir);
calc = calculate_equilibrium_voltage_from_soc(init_SOC);
[csn0, csp0] = init_cs(p, calc.Veq);
p.csn0 = csn0;
p.csp0 = csp0;

c_s_n0 = zeros(p.PadeOrder,1);
c_s_p0 = zeros(p.PadeOrder,1);

%%%%% Initial condition based on Jordan form
c_s_n0(3) = csn0;
c_s_p0(3) = csp0;

c_s_n = repmat(c_s_n0, [p.Nxn-1 1]);
c_s_p = repmat(c_s_p0, [p.Nxp-1 1]);

% Electrolyte concentration
c_e = p.c_e * ones(p.Nx-3,1);

% Spatial anode degradation states
delta_sei_n = zeros(p.Nxn-1,1);
delta_li_n = zeros(p.Nxn-1,1);
epsilon_e_n_dyn = p.epsilon_e_n * ones(p.Nxn-1,1);

% Temperature
T = t_amb + 273.15;

% Solid Potential
Uref_n0 = refPotentialAnode(p, csn0(1)*ones(p.Nxn-1,1) / p.c_s_n_max);
Uref_p0 = refPotentialCathode(p, csp0(1)*ones(p.Nxp-1,1) / p.c_s_p_max);

phi_s_n = Uref_n0;
phi_s_p = Uref_p0;

% Electrolyte Current
i_en = zeros(p.Nxn-1,1);
i_ep = zeros(p.Nxp-1,1);

% Electrolyte Potential
phi_e = zeros(p.Nx-1,1);

% Molar Ionic Flux
jn = zeros(p.Nxn-1,1);
jp = zeros(p.Nxp-1,1);

% Initial Conditions
if p.dfn_degradation_enabled
    p.x0 = [c_s_n; c_s_p; c_e; delta_sei_n; delta_li_n; epsilon_e_n_dyn; T];
else
    p.x0 = [c_s_n; c_s_p; c_e; T];
end
p.z0 = [phi_s_n; phi_s_p; i_en; i_ep; phi_e; jn; jp];

%%% Finite difference for spherical particle
p.Nr = Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;

% Finite difference points along x-coordinate
p.Nxn = Nxn;
p.Nxs = Nxs;
p.Nxp = Nxp;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

%% Output Discretization params
disp('Discretization Params');
fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
fprintf(1,'Order of Pade Approx for Solid Concentration : %1.0f\n',p.PadeOrder);
fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
disp(' ');

%% Set ambient temp and thermal dynamics switch (YK modified/added)
p.T_amb = t_amb + 273.15;
p.disable_thermal = disable_thermal;
p.delta_t = delta_t;

%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);

p.c_n_soc0 = cn_low;
p.c_p_soc0 = cp_low;
p.c_n_soc100 = cn_high;
p.c_p_soc100 = cp_high;
% Delta_cn = cn_high-cn_low;
% Delta_cp = cp_low-cp_high;
% p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
% p.OneC = p.OneC/p.Area_p; % Change to A/m^2 for DFN
%% Precompute data
% Vector lengths
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

% Solid concentration matrices
[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
p.A_csn = A_csn;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

clear A_csn B_csn A_csp B_csp C_csn C_csp;

% Electrolyte concentration matrices
[M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

rM3 = [Nn; Ns; Np];
cM3 = rM3';
p.ce.M3 = sparse(blkdiagFast(rM3, cM3, p.ce.M3n, p.ce.M3s, p.ce.M3p));

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce rM1 cM1;

% Solid Potential
[F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
    C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
p.F1_psn = F1_psn;
p.F1_psp = F1_psp;
p.F2_psn = F2_psn;
p.F2_psp = F2_psp;
p.G_psn = G_psn;
p.G_psp = G_psp;
p.C_psn = C_psn;
p.C_psp = C_psp;
p.D_psn = D_psn;
p.D_psp = D_psp;

clear F1_psn F1_psp F2_psn F2_psp G_psn G_psp C_psn C_psp D_psn D_psp;

% Electrolyte Current
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
p.F1_ien = F1_ien;
p.F1_iep = F1_iep;
p.F2_ien = F2_ien;
p.F2_iep = F2_iep;
p.F3_ien = F3_ien;
p.F3_iep = F3_iep;

clear F1_ien F1_iep F2_ien F2_iep F3_ien F3_iep;

% Electrolyte Potential
p.M1_pen_skel = sparse(diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
p.M1_pes_skel = sparse(diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
p.M1_pep_skel = sparse(diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

[M1_pe,M2_pe,M3_pe,M4_pe,C_pe] = phi_e_mats(p);
p.M1_pe = M1_pe;
p.M2_pe = M2_pe;
p.M3_pe = M3_pe;
p.M4_pe = M4_pe;
p.C_pe = C_pe;

clear M1_pe M2_pe M3_pe M4_pe C_pe

% Jacobian
[f_x, f_z, g_x, g_z] = jac_dfn_pre(p);
p.f_x = f_x;
p.f_z = f_z;
p.g_x = g_x;
p.g_z = g_z;
clear f_x f_z g_x g_z

% Initial Output States
[~, ~, p.y0] = dae_dfn(p.x0,p.z0,0,p);
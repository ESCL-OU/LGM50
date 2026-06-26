%% DAEs for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura

% Modified June 19, 2015 by Scott Moura
% Electrolyte concentration dynamics re-written with
%   - 2nd order boundary conditions
%   - rigorous implementation of nonlinear dynamics with Finite Difference

% Modified Apr 26, 2016 by Saehong Park
% Bruggeman relationship
%   - 'Effective' Electrolyte diffusion coefficients in here
%   - 'Effective' Electrolyte conductivity (Kappa) in phi_e_mats
%   - 'Effective' Solid-phase conductivity (Sigma) in phi_s_mats
%   - Modified electrolyte governing equation with volume fraction
%   - update phi_e_mats.m / phi_s_mats.m
%   - Added p.epsilon_f_n, p.epsilon_f_p in params_dualfoil


function [f, g, varargout] = dae_dfn(x,z,Cur,p)
import casadi.*
%% Compute A/m^2
I = Cur/p.Area_ref;
%% Parse out states
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
use_degradation = isfield(p, 'dfn_degradation_enabled') && p.dfn_degradation_enabled;
if use_degradation
    ind_delta_sei_n = Nc+1:Nc+Nn;
    ind_delta_li_n = Nc+Nn+1:Nc+2*Nn;
    ind_epsilon_e_n_dyn = Nc+2*Nn+1:Nc+3*Nn;
    ind_T = Nc+3*Nn+1;
else
    ind_T = Nc+1;
end
% x = dae_x(1:(Ncsn+Ncsp+Nx+1));
% z = dae_x((Ncsn+Ncsp+Nx+1+1):end);

% Solid Concentration
c_s_n = x(1:Ncsn);
c_s_p = x(Ncsn+1:Ncsn+Ncsp);

% Reformat into matrices
c_s_n_mat = reshape(c_s_n,p.PadeOrder,p.Nxn-1);
c_s_p_mat = reshape(c_s_p,p.PadeOrder,p.Nxp-1);

% Electrolyte concentration
c_e = x((Ncsn+Ncsp + 1):(Nc));

if use_degradation
    delta_sei_n = x(ind_delta_sei_n);
    delta_li_n = x(ind_delta_li_n);
    epsilon_e_n_dyn = x(ind_epsilon_e_n_dyn);
    epsilon_e_n = max(epsilon_e_n_dyn, p.epsilon_e_n_min);
else
    delta_sei_n = zeros(Nn,1);
    delta_li_n = zeros(Nn,1);
    epsilon_e_n = p.epsilon_e_n * ones(Nn,1);
end

% Temperature
if p.disable_thermal
    T = p.T_amb;
else
    T = x(ind_T);
end

% Solid Potential
phi_s_n = z(1:Nn);
phi_s_p = z(Nn+1:Nnp);

% Electrolyte Current
i_en = z(Nnp+1 : Nnp+Nn);
i_ep = z(Nnp+Nn+1 : 2*Nnp);

% Electrolyte Potential
phi_e = z(2*Nnp+1:2*Nnp+Nx+2);

% Molar ionic flux
jn = z(2*Nnp+Nx+3 : 2*Nnp+Nx+2+Nn);
jp = z(2*Nnp+Nx+2+Nn+1 : end);

%% Li Diffusion in Solid Phase: c_s(x,r,t)

% Anode LTI System
c_s_n_dot_mat = p.A_csn * c_s_n_mat + p.B_csn * jn';
y_csn = p.C_csn * c_s_n_mat;

% Anode Parse LTI outputs
c_ss_n = y_csn(1,:)';
c_avg_n = y_csn(2,:)';
c_s_n_dot = reshape(c_s_n_dot_mat,numel(c_s_n_dot_mat),1);

% Cathode LTI System
c_s_p_dot_mat = p.A_csp * c_s_p_mat + p.B_csp * jp';
y_csp = p.C_csp * c_s_p_mat;

% Cathode Parse LTI outputs
c_ss_p = y_csp(1,:)';
c_avg_p = y_csp(2,:)';
c_s_p_dot = reshape(c_s_p_dot_mat,numel(c_s_p_dot_mat),1);

%% Local anode side reactions
theta_n = c_ss_n / p.c_s_n_max;
theta_p = c_ss_p / p.c_s_p_max;
Unref = refPotentialAnode(p, theta_n);
Upref = refPotentialCathode(p, theta_p);

R_film_n = p.R_f_n + delta_sei_n ./ p.kappa_SEI;
phi_se_n = phi_s_n - phi_e(1:Nn);
eta_sei_n = phi_se_n - p.U_sei - p.Faraday * R_film_n .* jn;
if isfield(p, 'dfn_lpl_coupling_mode') && ~isempty(p.dfn_lpl_coupling_mode)
    lpl_coupling_mode = lower(char(p.dfn_lpl_coupling_mode));
else
    lpl_coupling_mode = 'spme_compatible';
end

j_sei_n = zeros(Nn,1);
delta_sei_dot_n = zeros(Nn,1);
if use_degradation
    sei_arg = (-p.alph*p.Faraday)/(p.R*T) .* eta_sei_n;
    sei_arg = min(max(sei_arg, p.dfn_exp_arg_min), p.dfn_exp_arg_max);
    K_sei = p.k0_sei .* exp(sei_arg);
    gamma_sei = p.M_sei/(2*p.rho_sei);
    delta_eff = max(delta_sei_n, 0) + p.delta_sei_min;
    A_bc = p.D_ec ./ delta_eff;
    B_bc = A_bc + K_sei;
    c_outer_ec = p.epsilon_sei * p.c0_ec;
    c_s_ec_n = (2*A_bc*c_outer_ec) ./ ...
        (B_bc + sqrt(B_bc.^2 + 4*gamma_sei.*K_sei.*A_bc*c_outer_ec));
    j_sei_n = -K_sei .* c_s_ec_n;
    delta_sei_dot_n = -gamma_sei .* j_sei_n;
end

j_lpl_n = zeros(Nn,1);
delta_li_dot_n = zeros(Nn,1);
film_flux_lpl_n = jn + j_sei_n;
if use_degradation && isfield(p, 'enable_lithium_plating') && p.enable_lithium_plating && Cur < 0
    % p.i0_lpl is an interfacial current density [A/m^2]; convert directly
    % to molar interfacial flux so M5n applies a_s only once.
    lpl_prefactor = -p.i0_lpl/p.Faraday;
    lpl_slope = (-p.alph_lpl*p.Faraday)/(p.R*T);
    switch lpl_coupling_mode
        case 'spme_compatible'
            % Default: match the existing SPMe-style side-reaction coupling.
            % Side-reaction drivers use the intercalation/SEI film drop while
            % electrolyte source/current balance still receives j_tot_n below.
            ldp_driver_n = phi_se_n - p.Faraday * R_film_n .* film_flux_lpl_n;
            lpl_arg = lpl_slope .* ldp_driver_n;
            lpl_arg = min(max(lpl_arg, p.dfn_exp_arg_min), p.dfn_exp_arg_max);
            j_lpl_n = lpl_prefactor .* exp(lpl_arg);
        case 'yang_total_flux'
            % Comparison mode: Yang-style implicit total reaction flux in
            % the LDP film drop, solved locally by fixed-point iteration.
            j_lpl_n = lpl_prefactor .* exp(min(max(lpl_slope .* ...
                (phi_se_n - p.Faraday * R_film_n .* film_flux_lpl_n), ...
                p.dfn_exp_arg_min), p.dfn_exp_arg_max));
            for kk = 1:p.dfn_lpl_max_iter
                j_old = j_lpl_n;
                ldp_driver_n = phi_se_n - p.Faraday * R_film_n .* (film_flux_lpl_n + j_old);
                lpl_arg = lpl_slope .* ldp_driver_n;
                lpl_arg = min(max(lpl_arg, p.dfn_exp_arg_min), p.dfn_exp_arg_max);
                j_lpl_n = lpl_prefactor .* exp(lpl_arg);
                if norm(j_lpl_n - j_old, inf) <= p.dfn_lpl_tol * max(1, norm(j_lpl_n, inf))
                    break;
                end
            end
        otherwise
            error('Unknown p.dfn_lpl_coupling_mode: %s.', lpl_coupling_mode);
    end
    delta_li_dot_n = -p.M_li/p.rho_li .* j_lpl_n;
    if any(abs(j_lpl_n) > p.dfn_lpl_warn_flux)
        warning('DFN:LiPlatingFluxLarge', ...
            'Large local Li plating flux detected: max(abs(j_lpl_n)) = %.3e mol/m^2-s.', ...
            max(abs(j_lpl_n)));
    end
end

j_tot_n = jn + j_sei_n + j_lpl_n;
if strcmp(lpl_coupling_mode, 'yang_total_flux')
    LDP_n = phi_se_n - p.Faraday * R_film_n .* j_tot_n;
else
    LDP_n = phi_se_n - p.Faraday * R_film_n .* film_flux_lpl_n;
end
LDP_total_n = phi_se_n - p.Faraday * R_film_n .* j_tot_n;
delta_film_dot_n = delta_sei_dot_n + delta_li_dot_n;
if use_degradation && isfield(p, 'enable_anode_porosity_loss') && p.enable_anode_porosity_loss
    epsilon_e_n_dot_raw = -p.a_s_n .* delta_film_dot_n;
    epsilon_e_n_dot = epsilon_e_n_dot_raw;
    epsilon_e_n_dot(epsilon_e_n_dyn <= p.epsilon_e_n_min & epsilon_e_n_dot_raw < 0) = 0;
else
    epsilon_e_n_dot = zeros(Nn,1);
end
if use_degradation && epsilon_e_n(end) <= 1.25*p.epsilon_e_n_min
    warning('DFN:AnodePorosityNearFloor', ...
        'Separator-side anode porosity is near floor: epsilon_e_n_dyn(end) = %.4g.', ...
        epsilon_e_n(end));
end

%% Li Diffusion in Electrolyte Phase: c_e(x,t)
% Separate and aggregate
c_en = c_e(1:Nn);
c_es = c_e(Nn+1:Nn+p.Nxs-1);
c_ep = c_e(Nn+p.Nxs : end);

% Compute electrolyte boundary concentrations. In degradation mode this
% rebuilds the anode/separator flux-continuity row online using the local
% separator-side anode porosity; this preserves the original matrix form
% while reflecting the dominant dynamic boundary effect.
if use_degradation
    Del_xn = p.L_n * p.delta_x_n;
    Del_xs = p.L_s * p.delta_x_s;
    Del_xp = p.L_p * p.delta_x_p;
    N1 = zeros(4,p.Nx-3);
    N2 = zeros(4);
    N1(1,1) = +4;
    N1(1,2) = -1;
    N2(1,1) = -3;
    eps_n_int = epsilon_e_n(end)^p.brug;
    eps_s = p.epsilon_e_s^p.brug;
    eps_p = p.epsilon_e_p^p.brug;
    N1(2,p.Nxn-2) = eps_n_int/(2*Del_xn);
    N1(2,p.Nxn-1) = -4*eps_n_int/(2*Del_xn);
    N2(2,2) = 3*eps_n_int/(2*Del_xn) + 3*eps_s/(2*Del_xs);
    N1(2,p.Nxn) = -4*eps_s/(2*Del_xs);
    N1(2,p.Nxn+1) = eps_s/(2*Del_xs);
    N1(3,p.Nxn+p.Nxs-3) = eps_s/(2*Del_xs);
    N1(3,p.Nxn+p.Nxs-2) = -4*eps_s/(2*Del_xs);
    N2(3,3) = 3*eps_s/(2*Del_xs) + 3*eps_p/(2*Del_xp);
    N1(3,p.Nxn+p.Nxs-1) = -4*eps_p/(2*Del_xp);
    N1(3,p.Nxn+p.Nxs) = eps_p/(2*Del_xp);
    N1(4,end-1) = +1;
    N1(4,end) = -4;
    N2(4,4) = +3;
    c_e_bcs = sparse(-N2\N1) * c_e;
else
    c_e_bcs = p.ce.C*c_e;
end
c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];
c_ex_log = max(c_ex, p.c_e_min);

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en,dD_en] = electrolyteDe(p, c_en, T);
[D_es,dD_es] = electrolyteDe(p, c_es, T);
[D_ep,dD_ep] = electrolyteDe(p, c_ep, T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
D_en_eff = D_en .* epsilon_e_n.^(p.brug-1);
dD_en_eff = dD_en .* epsilon_e_n.^(p.brug-1);

% DO same s,p 
D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO change D_en, D_es D_ep as D_eff... % Apr.22 2016 by Saehong Park

% System Matrices have all been precomputed & stored in param struct "p"

% Compute derivative
source_scale_n = p.epsilon_e_n ./ epsilon_e_n;
c_en_dot = dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) ...
    + p.ce.M5n*(source_scale_n .* j_tot_n);
if use_degradation && isfield(p, 'enable_anode_porosity_loss') && p.enable_anode_porosity_loss
    c_en_dot = c_en_dot - (epsilon_e_n_dot ./ epsilon_e_n) .* c_en;
end

c_es_dot = dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

c_ep_dot = dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + p.ce.M5p*jp;

% Assemble c_e_dot
c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

%% Potential in Solid Phase: phi_s(x,t)
% Algebraic eqns (semi-explicit form)
i_enn = [0; i_en; I];
i_epp = [I; i_ep; 0];

phi_sn_dot = p.F1_psn*phi_s_n + p.F2_psn*i_enn + p.G_psn*I;
phi_sp_dot = p.F1_psp*phi_s_p + p.F2_psp*i_epp + p.G_psp*I;

% Terminal Voltage
phi_s_n_bcs = p.C_psn * phi_s_n + p.D_psn * I;
phi_s_p_bcs = p.C_psp * phi_s_p + p.D_psp * I;

Volt = phi_s_p_bcs(2) - phi_s_n_bcs(1) - p.R_c*I;


%% Electrolyte Current: i_e(x,t)
i_en_dot = p.F1_ien*i_en + p.F2_ien*j_tot_n + p.F3_ien*I;
i_ep_dot = p.F1_iep*i_ep + p.F2_iep*jp + p.F3_iep*I;

% Electrolyte current across all three regions
i_e_in = [i_en; I*ones(p.Nxs+1,1); i_ep];

%% Potential in Electrolyte Phase: phi_e(x,t)

% Electrolyte Conductivity
kappa = electrolyteCond(p, c_ex_log, T);

kappa0 = kappa(1);                              % BC1
kappa_n = kappa(2:p.Nxn);
kappa_ns = kappa(p.Nxn+1);
kappa_s = kappa(p.Nxn+2 : p.Nxn+2+p.Nxs-2);
kappa_sp = kappa(p.Nxn+2+p.Nxs-1);
kappa_p = kappa(p.Nxn+2+p.Nxs : end-1);
kappaN = kappa(end);                            % BC2

% Effective conductivity - multiply by p.epsilon_e_x ^ (p.brug) % Apr.22 by Saehong Park
kappa_eff0 = kappa0 .* epsilon_e_n(1).^(p.brug);
kappa_eff_n = kappa_n .* epsilon_e_n.^(p.brug);

kappa_eff_ns = kappa_ns .* ((epsilon_e_n(end) + p.epsilon_e_s)/2).^(p.brug);

kappa_eff_s = kappa_s .* p.epsilon_e_s.^(p.brug);

kappa_eff_sp = kappa_sp .* ((p.epsilon_e_s + p.epsilon_e_p)/2).^(p.brug);

kappa_eff_p = kappa_p .* p.epsilon_e_p.^(p.brug);
kappa_effN = kappaN .* p.epsilon_e_p.^(p.brug);

% Form into vector
kappa_eff = [kappa_eff_n; kappa_eff_ns; kappa_eff_s; kappa_eff_sp; kappa_eff_p];
% Kap_eff = sparse(diag(kappa_eff));
Kap_eff = diag(kappa_eff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusional Conductivity
bet_org = (2*p.R*T)/(p.Faraday) * (p.t_plus - 1) * (1 + p.dactivity); 
Kap_eff_D = bet_org*Kap_eff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form Matrices
if isa(x, 'casadi.MX')
    M2_pe = MX(p.M2_pe);
else
    M2_pe = p.M2_pe;
end
M2_pe(1,1) = M2_pe(1,1) * kappa_eff0;
M2_pe(end,end) = M2_pe(end,end) * kappa_effN;

F1_pe = Kap_eff*p.M1_pe + M2_pe*p.C_pe;
F2_pe = p.M3_pe;
F3_pe = Kap_eff_D*p.M4_pe;

% Algebraic eqns (semi-explicit form)
phi_e_dot = F1_pe*phi_e + F2_pe*i_e_in + F3_pe*log(c_ex_log);

%% Butler-Volmer Equation
aFRT = (p.alph*p.Faraday)/(p.R*T);

% Exchange Current Density, i_0^{\pm}
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e,T);

% Overpotential, \eta
eta_n = phi_s_n - phi_e(1:Nn) - Unref - p.Faraday*R_film_n.*jn;
eta_p = phi_s_p - phi_e(end-Np+1:end) - Upref - p.Faraday*p.R_f_p*jp;

eta_pl = LDP_n;

% Algebraic eqns (semi-explicit form)
jn_dot = 2/p.Faraday * i_0n .* sinh(aFRT * eta_n) - jn;
jp_dot = 2/p.Faraday * i_0p .* sinh(aFRT * eta_p) - jp;

%% Temperature
if p.disable_thermal
    T_dot = 0;
else
    % Equilibrium Potential and Gradient wrt bulk concentration
    theta_avg_n = c_avg_n / p.c_s_n_max;
    theta_avg_p = c_avg_p / p.c_s_p_max;
    [Unb,~,dUnbdT] = refPotentialAnode(p, theta_avg_n);
    [Upb,~,dUpbdT] = refPotentialCathode(p, theta_avg_p);

    % Heat generated from intercalation (w/o boundaries for NOW)
    Q_nx = p.a_s_n*p.Faraday * jn .* (Unb - T*dUnbdT);
    Q_n = sum(Q_nx) * p.delta_x_n * p.L_n;

    Q_px = p.a_s_p*p.Faraday * jp .* (Upb - T*dUpbdT);
    Q_p = sum(Q_px) * p.delta_x_p * p.L_p;

    Q_inter = Q_n + Q_p;

    % Temperature ODE
    T_dot = (p.h*(p.T_amb - T) - I*Volt - Q_inter) / (p.rho_avg*p.C_p);
end

%% Concatenate Time Derivatives
if use_degradation
    f = [c_s_n_dot; c_s_p_dot; c_e_dot; ...
         delta_sei_dot_n; delta_li_dot_n; epsilon_e_n_dot; T_dot];
else
    f = [c_s_n_dot; c_s_p_dot; c_e_dot; T_dot];
end
g = [phi_sn_dot; phi_sp_dot; i_en_dot; i_ep_dot; ...
     phi_e_dot; jn_dot; jp_dot];

%% Optional Output Vars
% State of Charge
% TODO: Check SOC
SOC = (mean(c_avg_n) - p.c_n_soc0) / (p.c_n_soc100 - p.c_n_soc0);

% Conservation of Li-ion matters
% n_Li_s = sum(c_avg_n) * p.L_n*p.delta_x_n * p.epsilon_s_n * p.Area_n ...
%      + sum(c_avg_p) * p.L_p*p.delta_x_p * p.epsilon_s_p * p.Area_p;

% Total Li in solid phase
n_Li_s_n = p.Area_n * p.epsilon_s_n * sum(c_avg_n) * p.delta_x_n;
n_Li_s_p = p.Area_p * p.epsilon_s_p * sum(c_avg_p) * p.delta_x_p;
n_Li_s   = n_Li_s_n + n_Li_s_p;


n_Li_e = sum(c_e(1:Nn) .* epsilon_e_n) * p.L_n*p.delta_x_n * p.Area_n ...
     + sum(c_e(Nn+1:end-Np)) * p.L_s*p.delta_x_s * p.epsilon_e_s * p.Area_s ...
     + sum(c_e(end-Np+1:end)) * p.L_p*p.delta_x_p * p.epsilon_e_p * p.Area_p;

% nLidot = sum(jn) * p.delta_x_n * Nn + sum(jp) * p.delta_x_p * Np;

c_e0n = c_ex(1);
c_e0p = c_ex(end);

eta_s_Ln = phi_s_n_bcs(2) - phi_e(Nn+1);

ce0n = c_e_bcs(1);
cens = c_e_bcs(2);
cesp = c_e_bcs(3);
ce0p = c_e_bcs(4);

% Separate and aggregate electrolyte concentration
c_en = c_e(1:(p.Nxn-1));
c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
c_ep = c_e((p.Nxn-1)+p.Nxs : end);
c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];

% Aggregate Outputs
dfn_deg_y = [LDP_n; LDP_total_n; jn; j_sei_n; j_lpl_n; j_tot_n; ...
             delta_sei_n; delta_li_n; epsilon_e_n; phi_e(1:Nn); phi_s_n; ...
             LDP_n(end); LDP_total_n(end); j_lpl_n(end); delta_li_n(end); epsilon_e_n(end)];
if use_degradation
    varargout{1} = [Volt; SOC; eta_pl; c_ss_n; c_ss_p; c_avg_n; c_avg_p; c_ex; eta_n; eta_p; Unref; dfn_deg_y];
else
    varargout{1} = [Volt; SOC; eta_pl; c_ss_n; c_ss_p; c_avg_n; c_avg_p; c_ex; eta_n; eta_p; Unref];
end

if nargout > 3
    IR_n = p.Faraday*R_film_n.*jn;
    IR_p = p.Faraday*p.R_f_p*jp;
    phi_e_n = phi_e(1:Nn);
    phi_e_p = phi_e(end-Np+1:end);

%     V_noVCE = Volt - (phi_e(end) - phi_e(1));
    V_noVCE = eta_p(end) - eta_n(1) ...
            + Upref(end) - Unref(1) ...
            - (IR_p(end) + IR_n(1));
    
    % Overpotential due to electrolyte conductivity
    V_electrolyteCond = (p.L_n/(2*mean(kappa_eff_n)) + 2*p.L_s/(2*mean(kappa_eff_s)) + p.L_p/(2*mean(kappa_eff_p)))*Cur;

    % Overpotential due to electrolyte polarization
    V_electrolytePolar = (2*p.R*p.T_amb)/(p.Faraday) * ...
                            (1 * (log(max(cens,p.c_e_min)) - log(max(ce0n,p.c_e_min))) ...
                            +1 * (log(max(cesp,p.c_e_min)) - log(max(cens,p.c_e_min))) ...
                            +1 * (log(max(ce0p,p.c_e_min)) - log(max(cesp,p.c_e_min))));

    V_electrolyte = phi_e_p(end) - phi_e_n(1);

    varargout{2} = eta_n;
    varargout{3} = eta_p;
    varargout{4} = Unref;
    varargout{5} = Upref;
    varargout{6} = IR_n;
    varargout{7} = IR_p;
    varargout{8} = phi_e_n;
    varargout{9} = phi_e_p;
    varargout{10} = V_noVCE;
    varargout{11} = V_electrolyteCond;
    varargout{12} = V_electrolytePolar;
    varargout{13} = LDP_n;
    varargout{14} = LDP_total_n;
    varargout{15} = jn;
    varargout{16} = j_sei_n;
    varargout{17} = j_lpl_n;
    varargout{18} = j_tot_n;
    varargout{19} = delta_sei_n;
    varargout{20} = delta_li_n;
    varargout{21} = epsilon_e_n;
    varargout{22} = phi_s_n;
end


function out = test_spme_legacy_vs_current()
%TEST_SPME_LEGACY_VS_CURRENT
% Compare voltage response of:
%   - current SPMe.m
%   - recreated legacy SPMe_legacy.m (pre-update)
% with SEI disabled.

    tests_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(tests_dir);
    addpath(repo_dir);
    addpath(fullfile(repo_dir, 'SPMe'));
    addpath(fullfile(repo_dir, 'param'));

    % Test condition (CC charge from V_min to V_max).
    init_SOC = 0.0;
    CRate = -1.0;
    t_amb = 25;
    delta_t = 1;
    Nr = 20; Nxn = 20; Nxs = 10; Nxp = 20;
    disable_thermal = true;

    sim_new = sim_profile_LGM50('SPMe', init_SOC, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'none', 'SaveOutput', false);
    p = sim_new.p;
    I = sim_new.current(1);
    t = sim_new.time(:);
    NT = numel(t);

    m_old = SPMe_legacy(p, 'idas');
    x_old = p.x0;
    y0_old = m_old.h('x', x_old, 'p', I); y0_old = full(y0_old.y);

    v_new = sim_new.voltage(:);
    v_old = zeros(NT,1);
    v_noVCE_new = zeros(NT,1);
    v_noVCE_old = zeros(NT,1);
    v_electrolyteCond_new = zeros(NT,1);
    v_electrolyteCond_old = zeros(NT,1);
    v_electrolytePolar_new = zeros(NT,1);
    v_electrolytePolar_old = zeros(NT,1);

    v_old(1) = y0_old(1);
    [v_noVCE_new(1), v_electrolyteCond_new(1), v_electrolytePolar_new(1)] = calc_spme_voltage_terms(p, sim_new.x(:,1), I);
    [~, ~, ~, ~, v_noVCE_old(1), v_electrolyteCond_old(1), v_electrolytePolar_old(1)] = m_old.potentials(x_old, I);

    t_last = NT;
    for k = 1:NT-1
        [x_old, y_old] = m_old.step(x_old, I);
        v_old(k+1) = y_old(1);
        [v_noVCE_new(k+1), v_electrolyteCond_new(k+1), v_electrolytePolar_new(k+1)] = calc_spme_voltage_terms(p, sim_new.x(:,k+1), I);
        [~, ~, ~, ~, v_noVCE_old(k+1), v_electrolyteCond_old(k+1), v_electrolytePolar_old(k+1)] = m_old.potentials(x_old, I);
    end

    t = t(1:t_last);
    v_new = v_new(1:t_last);
    v_old = v_old(1:t_last);
    v_noVCE_new = v_noVCE_new(1:t_last);
    v_noVCE_old = v_noVCE_old(1:t_last);
    v_electrolyteCond_new = v_electrolyteCond_new(1:t_last);
    v_electrolyteCond_old = v_electrolyteCond_old(1:t_last);
    v_electrolytePolar_new = v_electrolytePolar_new(1:t_last);
    v_electrolytePolar_old = v_electrolytePolar_old(1:t_last);
    dv = v_new - v_old;
    rmse_v = sqrt(mean(dv.^2));

    figure('Name', 'SPMe current vs legacy (SEI disabled)', 'Color', 'w');
    subplot(5,1,1);
    plot(t, v_new, 'b-', 'LineWidth', 1.4); hold on;
    plot(t, v_old, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    title(sprintf('Voltage comparison (SEI disabled), RMSE = %.6g V', rmse_v));
    legend('Current SPMe', 'Legacy SPMe', 'Location', 'best');

    subplot(5,1,2);
    plot(t, dv, 'k-', 'LineWidth', 1.2); grid on;
    ylabel('\DeltaV = V_{new} - V_{legacy} [V]');
    title('Voltage difference');

    subplot(5,1,3);
    plot(t, v_noVCE_new, 'b-', 'LineWidth', 1.4); hold on;
    plot(t, v_noVCE_old, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('V_{noVCE} [V]');
    title('No-VCE term');

    subplot(5,1,4);
    plot(t, v_electrolyteCond_new, 'b-', 'LineWidth', 1.4); hold on;
    plot(t, v_electrolyteCond_old, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('V_{elyteCond} [V]');
    title('Electrolyte conductivity term');

    subplot(5,1,5);
    plot(t, v_electrolytePolar_new, 'b-', 'LineWidth', 1.4); hold on;
    plot(t, v_electrolytePolar_old, 'r--', 'LineWidth', 1.2); grid on;
    xlabel('Time [s]');
    ylabel('V_{elytePolar} [V]');
    title('Electrolyte polarization term');

    out = struct();
    out.time = t;
    out.v_new = v_new;
    out.v_legacy = v_old;
    out.dv = dv;
    out.v_noVCE_new = v_noVCE_new;
    out.v_noVCE_legacy = v_noVCE_old;
    out.v_electrolyteCond_new = v_electrolyteCond_new;
    out.v_electrolyteCond_legacy = v_electrolyteCond_old;
    out.v_electrolytePolar_new = v_electrolytePolar_new;
    out.v_electrolytePolar_legacy = v_electrolytePolar_old;
    out.rmse_v = rmse_v;
    out.conditions = struct('init_SOC', init_SOC, 'CRate', CRate, 't_amb', t_amb, ...
                            'delta_t', delta_t, 'degradation_model', p.degradation_model);
end

function [V_noVCE, V_electrolyteCond, V_electrolytePolar] = calc_spme_voltage_terms(p, x, I)
    c_s_n = x(1:p.Nr-1);
    c_s_p = x(p.Nr : 2*(p.Nr-1));
    c_e = x(2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);

    jn_tot = I/(p.Faraday*p.a_s_n*p.Area_n*p.L_n);
    jp_tot = -I/(p.Faraday*p.a_s_p*p.Area_p*p.L_p);
    [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_obs_mats(p, p.T_amb); %#ok<ASGLU>
    c_ss_n = C_n*c_s_n + D_n*jn_tot;
    c_ss_p = C_p*c_s_p + D_p*jp_tot;

    [~, c_ex] = ode_electrolyte(c_e, jn_tot, jp_tot, p);
    ce0n = c_ex(1);
    cens = c_ex(p.Nxn+1);
    cesp = c_ex(p.Nxn+p.Nxs+1);
    ce0p = c_ex(end);
    cen_bar = mean(c_ex(1:p.Nxn+1));
    ces_bar = mean(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1)));
    cep_bar = mean(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1)));

    ce_bar = zeros(p.Nx, 1);
    ce_bar(1:p.Nxn) = cen_bar;
    ce_bar((p.Nxn+1):(p.Nxn+p.Nxs-1)) = ces_bar;
    ce_bar(p.Nxn+p.Nxs+1:end) = cep_bar;
    [i_0n,i_0p] = exch_cur_dens(p, c_ss_n, c_ss_p, ce_bar, p.T_amb);
    RTaF = (p.R*p.T_amb)/(p.alph*p.Faraday);
    Unref = refPotentialAnode(p, c_ss_n/p.c_s_n_max);
    Upref = refPotentialCathode(p, c_ss_p/p.c_s_p_max);
    eta_n = RTaF * asinh(I / (2*p.a_s_n*p.Area_n*p.L_n*i_0n(1)));
    eta_p = RTaF * asinh(-I / (2*p.a_s_p*p.Area_p*p.L_p*i_0p(end)));

    V_noVCE = eta_p - eta_n + Upref - Unref ...
        - (p.R_f_n/(p.a_s_n*p.L_n*p.Area_n) + p.R_f_p/(p.a_s_p*p.L_p*p.Area_p))*I;

    kap_n = electrolyteCond(p, cen_bar, p.T_amb);
    kap_s = electrolyteCond(p, ces_bar, p.T_amb);
    kap_p = electrolyteCond(p, cep_bar, p.T_amb);
    kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
    kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
    kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);
    V_electrolyteCond = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*I;

    dfca_n = electrolyteAct(cen_bar, p.T_amb, p);
    dfca_s = electrolyteAct(ces_bar, p.T_amb, p);
    dfca_p = electrolyteAct(cep_bar, p.T_amb, p);
    V_electrolytePolar = (2*p.R*p.T_amb)/(p.Faraday) * ...
        ((1+dfca_n) * (log(cens) - log(ce0n)) + ...
         (1+dfca_s) * (log(cesp) - log(cens)) + ...
         (1+dfca_p) * (log(ce0p) - log(cesp)));
end

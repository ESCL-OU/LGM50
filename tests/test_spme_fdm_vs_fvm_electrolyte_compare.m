function out = test_spme_fdm_vs_fvm_electrolyte_compare()
%TEST_SPME_FDM_VS_FVM_ELECTROLYTE_COMPARE
% Compare SPMe electrolyte dynamics implementations:
%   - fdm (existing discretization)
%   - fvm (conservative discretization)
%
% Uses the same sequence profile style as
% test_sim_profile_sequence_capacity_cycle_plot and focuses on the same
% c_e-related diagnostics:
%   1) per-step lithium inventory change terms
%   2) electrolyte inventory by region and total

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    % ----------------------------
    % Adjustable discretization
    % ----------------------------
    Nxn = 10;
    Nxs = 5;
    Nxp = 10;
    Nr = 15;
    delta_t = 1;
    t_amb = 25;

    % ----------------------------
    % Test profile (same style as your sequence test)
    % ----------------------------
    model_name = 'SPMe';
    init_soc = 0.0;
    sequence = {'REST', 0.5, 'CC', -1, 'REST', 0.5, 'CC', 1};
    repeat_count = 1;

    req_outputs = {'cen_bar', 'ces_bar', 'cep_bar', ...
        'cs_ave_n_dot', 'cs_ave_p_dot', ...
        'c_en_dot_bar', 'c_es_dot_bar', 'c_ep_dot_bar'};

    sim_fdm = sim_profile_sequence_LGM50(model_name, init_soc, sequence, ...
        'Repeat', repeat_count, ...
        'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, 'Nr', Nr, ...
        'DeltaT', delta_t, 'Tamb', t_amb, ...
        'DegradationModel', 'none', ...
        'ElectrolyteModel', 'fdm', ...
        'SPMeOutputs', req_outputs, ...
        'SaveOutput', false, 'ShowProgress', false);

    sim_fvm = sim_profile_sequence_LGM50(model_name, init_soc, sequence, ...
        'Repeat', repeat_count, ...
        'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, 'Nr', Nr, ...
        'DeltaT', delta_t, 'Tamb', t_amb, ...
        'DegradationModel', 'none', ...
        'ElectrolyteModel', 'fvm', ...
        'SPMeOutputs', req_outputs, ...
        'SaveOutput', false, 'ShowProgress', false);

    m_fdm = calc_ce_metrics(sim_fdm);
    m_fvm = calc_ce_metrics(sim_fvm);

    figure('Name', 'SPMe electrolyte compare: FDM vs FVM', 'Color', 'w');
    tl = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(1);
    plot(m_fdm.time, m_fdm.dN_e_n_step, 'b-', 'LineWidth', 0.9); hold on;
    plot(m_fdm.time, m_fdm.dN_e_s_step, 'm-', 'LineWidth', 0.9);
    plot(m_fdm.time, m_fdm.dN_e_p_step, 'r-', 'LineWidth', 0.9);
    plot(m_fdm.time, m_fdm.dN_s_n_step, 'c-', 'LineWidth', 0.9);
    plot(m_fdm.time, m_fdm.dN_s_p_step, 'y-', 'LineWidth', 0.9);
    plot(m_fdm.time, m_fdm.dN_total_step, 'k-', 'LineWidth', 1.2);
    grid on;
    ylabel('\DeltaN_{Li} per step [mol]');
    title('FDM: per-step lithium inventory terms');
    legend('\DeltaN_{e,n}', '\DeltaN_{e,s}', '\DeltaN_{e,p}', ...
        '\DeltaN_{s,n}', '\DeltaN_{s,p}', '\DeltaN_{total}', ...
        'Location', 'best');

    nexttile(2);
    plot(m_fvm.time, m_fvm.dN_e_n_step, 'b-', 'LineWidth', 0.9); hold on;
    plot(m_fvm.time, m_fvm.dN_e_s_step, 'm-', 'LineWidth', 0.9);
    plot(m_fvm.time, m_fvm.dN_e_p_step, 'r-', 'LineWidth', 0.9);
    plot(m_fvm.time, m_fvm.dN_s_n_step, 'c-', 'LineWidth', 0.9);
    plot(m_fvm.time, m_fvm.dN_s_p_step, 'y-', 'LineWidth', 0.9);
    plot(m_fvm.time, m_fvm.dN_total_step, 'k-', 'LineWidth', 1.2);
    grid on;
    ylabel('\DeltaN_{Li} per step [mol]');
    title('FVM: per-step lithium inventory terms');
    legend('\DeltaN_{e,n}', '\DeltaN_{e,s}', '\DeltaN_{e,p}', ...
        '\DeltaN_{s,n}', '\DeltaN_{s,p}', '\DeltaN_{total}', ...
        'Location', 'best');

    nexttile(3);
    plot(m_fdm.time, m_fdm.n_e_n, 'b-', 'LineWidth', 1.0); hold on;
    plot(m_fdm.time, m_fdm.n_e_s, 'm-', 'LineWidth', 1.0);
    plot(m_fdm.time, m_fdm.n_e_p, 'r-', 'LineWidth', 1.0);
    plot(m_fdm.time, m_fdm.n_e_total, 'k-', 'LineWidth', 1.3);
    grid on;
    ylabel('N_e [mol]');
    xlabel('Time [s]');
    title('FDM: electrolyte inventory');
    legend('N_{e,n}', 'N_{e,s}', 'N_{e,p}', 'N_{e,total}', ...
        'Location', 'best');

    nexttile(4);
    plot(m_fvm.time, m_fvm.n_e_n, 'b-', 'LineWidth', 1.0); hold on;
    plot(m_fvm.time, m_fvm.n_e_s, 'm-', 'LineWidth', 1.0);
    plot(m_fvm.time, m_fvm.n_e_p, 'r-', 'LineWidth', 1.0);
    plot(m_fvm.time, m_fvm.n_e_total, 'k-', 'LineWidth', 1.3);
    grid on;
    ylabel('N_e [mol]');
    xlabel('Time [s]');
    title('FVM: electrolyte inventory');
    legend('N_{e,n}', 'N_{e,s}', 'N_{e,p}', 'N_{e,total}', ...
        'Location', 'best');

    title(tl, sprintf('SPMe FDM vs FVM (Nxn=%d, Nxs=%d, Nxp=%d, Nr=%d)', ...
        Nxn, Nxs, Nxp, Nr));

    out = struct();
    out.settings = struct('Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, 'Nr', Nr, ...
        'DeltaT', delta_t, 'Tamb_degC', t_amb, 'init_soc', init_soc, ...
        'sequence', {sequence}, 'repeat', repeat_count);
    out.fdm = struct('sim', sim_fdm, 'metrics', m_fdm);
    out.fvm = struct('sim', sim_fvm, 'metrics', m_fvm);
end

function m = calc_ce_metrics(sim)
    p = sim.p;
    dt = p.delta_t;

    cen_bar = sim.spme_outputs.cen_bar(:);
    ces_bar = sim.spme_outputs.ces_bar(:);
    cep_bar = sim.spme_outputs.cep_bar(:);
    c_en_dot_bar = sim.spme_outputs.c_en_dot_bar(:);
    c_es_dot_bar = sim.spme_outputs.c_es_dot_bar(:);
    c_ep_dot_bar = sim.spme_outputs.c_ep_dot_bar(:);
    cs_ave_n_dot = sim.spme_outputs.cs_ave_n_dot(:);
    cs_ave_p_dot = sim.spme_outputs.cs_ave_p_dot(:);

    vol_e_n = p.epsilon_e_n * p.Area_n * p.L_n;
    vol_e_s = p.epsilon_e_s * p.Area_s * p.L_s;
    vol_e_p = p.epsilon_e_p * p.Area_p * p.L_p;
    vol_s_n = p.epsilon_s_n * p.Area_n * p.L_n;
    vol_s_p = p.epsilon_s_p * p.Area_p * p.L_p;

    m = struct();
    m.time = sim.time(:);
    m.n_e_n = cen_bar * vol_e_n;
    m.n_e_s = ces_bar * vol_e_s;
    m.n_e_p = cep_bar * vol_e_p;
    m.n_e_total = m.n_e_n + m.n_e_s + m.n_e_p;

    m.dN_e_n_step = c_en_dot_bar * vol_e_n * dt;
    m.dN_e_s_step = c_es_dot_bar * vol_e_s * dt;
    m.dN_e_p_step = c_ep_dot_bar * vol_e_p * dt;
    m.dN_s_n_step = cs_ave_n_dot * vol_s_n * dt;
    m.dN_s_p_step = cs_ave_p_dot * vol_s_p * dt;
    m.dN_total_step = m.dN_e_n_step + m.dN_e_s_step + m.dN_e_p_step ...
        + m.dN_s_n_step + m.dN_s_p_step;
end

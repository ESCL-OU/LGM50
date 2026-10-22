function out = test_capacity_LGM50()
%TEST_CAPACITY_LGM50 Capacity check from params_LGM50 and SPMe discharge.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    param_dir = fullfile(base_dir, 'param');
    spme_dir = fullfile(base_dir, 'SPMe');

    if ~isfolder(param_dir) || ~isfolder(spme_dir)
        error('Required folders ''param'' and ''SPMe'' were not found under %s.', base_dir);
    end

    addpath(base_dir);
    addpath(param_dir);
    addpath(spme_dir);

    run(fullfile(param_dir, 'params_LGM50.m'));

    [cn_low, cp_low] = init_cs(p, p.volt_min);
    [cn_high, cp_high] = init_cs(p, p.volt_max);

    nLi_anode = (cn_high - cn_low) * p.epsilon_s_n * p.L_n * p.Area_n;
    nLi_cathode = (cp_low - cp_high) * p.epsilon_s_p * p.L_p * p.Area_p;

    Q_anode_Ah = nLi_anode * p.Faraday / 3600;
    Q_cathode_Ah = nLi_cathode * p.Faraday / 3600;
    Q_theory_Ah = 0.5 * (Q_anode_Ah + Q_cathode_Ah);

    d_nLi = nLi_anode - nLi_cathode;
    li_balance_rel_pct = 100 * d_nLi / max(abs(nLi_anode), eps);

    CRate = 1.0;
    t_amb = 25;
    delta_t = 1;
    Nr = 20;
    Nxn = 20;
    Nxs = 10;
    Nxp = 20;

    sim_spme = sim_profile_LGM50('SPMe', 1.0, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'none', ...
        'ElectrolyteModel', 'fvm', ...
        'SPMeOutputs', {'cs_ave_n', 'cs_ave_p', 'cen_bar', 'ces_bar', 'cep_bar'}, ...
        'SaveOutput', false);
    p_spme = sim_spme.p;
    y0 = sim_spme.voltage(1);
    t_cutoff = sim_spme.time(end);
    Q_sim_Ah = trapz(sim_spme.time, sim_spme.current) / 3600;
    sim_vs_theory_rel_pct = 100 * (Q_sim_Ah - Q_theory_Ah) / max(Q_theory_Ah, eps);
    I_cmd = p_spme.OneC * CRate;
    I_mean = mean(sim_spme.current);
    I_min = min(sim_spme.current);
    I_max = max(sim_spme.current);
    I_std = std(sim_spme.current);
    C_rate_mean = I_mean / max(p_spme.OneC, eps);
    C_rate_min = I_min / max(p_spme.OneC, eps);
    C_rate_max = I_max / max(p_spme.OneC, eps);
    Q_from_mean_Ah = I_mean * t_cutoff / 3600;
    Ah_cum = cumtrapz(sim_spme.time, sim_spme.current) / 3600;

    cs_ave_n = sim_spme.spme_outputs.cs_ave_n(:);
    cs_ave_p = sim_spme.spme_outputs.cs_ave_p(:);
    cen_bar = sim_spme.spme_outputs.cen_bar(:);
    ces_bar = sim_spme.spme_outputs.ces_bar(:);
    cep_bar = sim_spme.spme_outputs.cep_bar(:);

    vol_s_n = p_spme.epsilon_s_n * p_spme.Area_n * p_spme.L_n;
    vol_s_p = p_spme.epsilon_s_p * p_spme.Area_p * p_spme.L_p;
    vol_e_n = p_spme.epsilon_e_n * p_spme.Area_n * p_spme.L_n;
    vol_e_s = p_spme.epsilon_e_s * p_spme.Area_s * p_spme.L_s;
    vol_e_p = p_spme.epsilon_e_p * p_spme.Area_p * p_spme.L_p;

    k0 = 1;
    kf = numel(sim_spme.time);
    dNs_n_mol = (cs_ave_n(kf) - cs_ave_n(k0)) * vol_s_n;
    dNs_p_mol = (cs_ave_p(kf) - cs_ave_p(k0)) * vol_s_p;
    dNe_mol = (cen_bar(kf) - cen_bar(k0)) * vol_e_n + ...
        (ces_bar(kf) - ces_bar(k0)) * vol_e_s + ...
        (cep_bar(kf) - cep_bar(k0)) * vol_e_p;
    li_residual_mol = dNs_n_mol + dNs_p_mol + dNe_mol;
    li_residual_Ah = li_residual_mol * p_spme.Faraday / 3600;

    % For discharge profile: anode decreases, cathode increases.
    Q_profile_anode_Ah = max(cs_ave_n(k0) - cs_ave_n(kf), 0) * vol_s_n * p_spme.Faraday / 3600;
    Q_profile_cathode_Ah = max(cs_ave_p(kf) - cs_ave_p(k0), 0) * vol_s_p * p_spme.Faraday / 3600;
    sim_vs_profile_anode_rel_pct = 100 * (Q_sim_Ah - Q_profile_anode_Ah) / max(Q_profile_anode_Ah, eps);
    sim_vs_profile_cathode_rel_pct = 100 * (Q_sim_Ah - Q_profile_cathode_Ah) / max(Q_profile_cathode_Ah, eps);
    anode_vs_cathode_Ah = Q_profile_anode_Ah - Q_profile_cathode_Ah;
    anode_vs_cathode_rel_pct = 100 * anode_vs_cathode_Ah / max(Q_profile_anode_Ah, eps);

    fprintf("\n=== Capacity from init_cs + params_LGM50 ===\n");
    fprintf("cn(Vmin): %.4f | cn(Vmax): %.4f [mol/m^3]\n", cn_low, cn_high);
    fprintf("cp(Vmin): %.4f | cp(Vmax): %.4f [mol/m^3]\n", cp_low, cp_high);
    fprintf("Moved Li (anode decrease)  : %.8f mol\n", nLi_anode);
    fprintf("Moved Li (cathode increase): %.8f mol\n", nLi_cathode);
    fprintf("Li balance mismatch (an-cath): %.3e mol (%.4f %%)\n", d_nLi, li_balance_rel_pct);
    fprintf("Capacity from anode   : %.6f Ah\n", Q_anode_Ah);
    fprintf("Capacity from cathode : %.6f Ah\n", Q_cathode_Ah);
    fprintf("Theoretical capacity (avg): %.6f Ah\n", Q_theory_Ah);

    fprintf("\n=== SPMe discharge check (%.3gC) ===\n", CRate);
    fprintf("Initial OCV at SPMe start: %.6f V (target Vmax %.6f V)\n", y0, p_spme.volt_max);
    fprintf("Cutoff time to Vmin       : %.2f s\n", t_cutoff);
    fprintf("OneC current (p.OneC)     : %.6f A\n", p_spme.OneC);
    fprintf("Commanded current (%.3gC)  : %.6f A\n", CRate, I_cmd);
    fprintf("Measured current stats     : mean=%.6f A, min=%.6f A, max=%.6f A, std=%.3e A\n", ...
        I_mean, I_min, I_max, I_std);
    fprintf("Measured C-rate stats      : mean=%.6fC, min=%.6fC, max=%.6fC\n", ...
        C_rate_mean, C_rate_min, C_rate_max);
    fprintf("Capacity from mean I*time  : %.6f Ah\n", Q_from_mean_Ah);
    fprintf("Delivered capacity (SPMe) : %.6f Ah\n", Q_sim_Ah);
    fprintf("SPMe vs theory error      : %.4f %%\n", sim_vs_theory_rel_pct);

    fprintf("\n=== SPMe profile-window capacity + Li-balance diagnostics ===\n");
    fprintf("Window: start of profile -> end of profile (k=%d -> %d)\n", k0, kf);
    fprintf("Q_profile_anode (from cs_ave_n)   : %.6f Ah\n", Q_profile_anode_Ah);
    fprintf("Q_profile_cathode (from cs_ave_p) : %.6f Ah\n", Q_profile_cathode_Ah);
    fprintf("Q_integral (trapz current)        : %.6f Ah\n", Q_sim_Ah);
    fprintf("Integral vs anode error           : %.4f %%\n", sim_vs_profile_anode_rel_pct);
    fprintf("Integral vs cathode error         : %.4f %%\n", sim_vs_profile_cathode_rel_pct);
    fprintf("Anode-cathode difference          : %+.6e Ah (%.4f %%)\n", ...
        anode_vs_cathode_Ah, anode_vs_cathode_rel_pct);
    fprintf("dNs_n=%+.6e mol, dNs_p=%+.6e mol, dNe=%+.6e mol\n", ...
        dNs_n_mol, dNs_p_mol, dNe_mol);
    fprintf("Residual dNs_n+dNs_p+dNe          : %+.6e mol (%+.6e Ah)\n\n", ...
        li_residual_mol, li_residual_Ah);

    figure('Name', 'test\_capacity\_LGM50 diagnostics', 'Color', 'w');
    tl = tiledlayout(4,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(sim_spme.time, sim_spme.current, 'm-', 'LineWidth', 1.1); hold on;
    yline(I_cmd, '--k', 'I_{cmd}', 'LineWidth', 0.9);
    grid on;
    ylabel('Current [A]');
    legend('I', 'I_{cmd}', 'Location', 'best');
    title(sprintf('SPMe CC discharge: commanded %.3gC, measured mean %.4fC', ...
        CRate, C_rate_mean));

    nexttile;
    plot(sim_spme.time, sim_spme.voltage, 'k-', 'LineWidth', 1.1);
    grid on;
    ylabel('Voltage [V]');

    nexttile;
    plot(sim_spme.time, sim_spme.soc, 'b-', 'LineWidth', 1.1);
    grid on;
    ylabel('SoC [-]');

    nexttile;
    plot(sim_spme.time, Ah_cum, 'r-', 'LineWidth', 1.2); hold on;
    yline(Q_sim_Ah, '--k', sprintf('Q=%.4f Ah', Q_sim_Ah), 'LineWidth', 0.9);
    grid on;
    xlabel('Time [s]');
    ylabel('Cum. Ah [Ah]');
    title(tl, 'test\_capacity\_LGM50');

    out = struct();
    out.cn_low = cn_low;
    out.cn_high = cn_high;
    out.cp_low = cp_low;
    out.cp_high = cp_high;
    out.nLi_anode = nLi_anode;
    out.nLi_cathode = nLi_cathode;
    out.li_balance_mismatch_mol = d_nLi;
    out.li_balance_mismatch_rel_pct = li_balance_rel_pct;
    out.Q_anode_Ah = Q_anode_Ah;
    out.Q_cathode_Ah = Q_cathode_Ah;
    out.Q_theory_Ah = Q_theory_Ah;
    out.Q_spme_Ah = Q_sim_Ah;
    out.Q_from_mean_Ah = Q_from_mean_Ah;
    out.spme_vs_theory_rel_pct = sim_vs_theory_rel_pct;
    out.Q_profile_anode_Ah = Q_profile_anode_Ah;
    out.Q_profile_cathode_Ah = Q_profile_cathode_Ah;
    out.spme_vs_profile_anode_rel_pct = sim_vs_profile_anode_rel_pct;
    out.spme_vs_profile_cathode_rel_pct = sim_vs_profile_cathode_rel_pct;
    out.anode_vs_cathode_Ah = anode_vs_cathode_Ah;
    out.anode_vs_cathode_rel_pct = anode_vs_cathode_rel_pct;
    out.dNs_n_mol = dNs_n_mol;
    out.dNs_p_mol = dNs_p_mol;
    out.dNe_mol = dNe_mol;
    out.li_residual_mol = li_residual_mol;
    out.li_residual_Ah = li_residual_Ah;
    out.current_command_A = I_cmd;
    out.current_mean_A = I_mean;
    out.current_min_A = I_min;
    out.current_max_A = I_max;
    out.current_std_A = I_std;
    out.c_rate_mean = C_rate_mean;
    out.c_rate_min = C_rate_min;
    out.c_rate_max = C_rate_max;
    out.cumulative_Ah = Ah_cum(:);
    out.spme_initial_voltage = y0;
    out.spme_cutoff_time_s = t_cutoff;
end

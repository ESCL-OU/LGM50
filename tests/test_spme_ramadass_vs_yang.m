function out = test_spme_ramadass_vs_yang()
%TEST_SPME_RAMADASS_VS_YANG Compare SPMe degradation quantities.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));
    results_dir = fullfile(base_dir, 'results');
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    run_ts = datestr(now, 'yyyymmdd_HHMMSS');

    t_amb = 25;
    delta_t = 1;
    init_soc_charge = 0.0;
    sequence1 = {'REST', 0.5, 'CC', -0.2, 'REST', 0.5, 'CC', 0.2};
    sequence2 = {'CC', -1.0, 'CC', 1.0};
    sequence1_repeat = 1;
    sequence2_repeat = 3;
    block_sequence = [repmat(sequence1, 1, sequence1_repeat), ...
        repmat(sequence2, 1, sequence2_repeat)];
    repeat_count = 2;
    sequence = repmat(block_sequence, 1, repeat_count);
    Nr = 15; Nxn = 30; Nxs = 5; Nxp = 30;

    logmsg('Starting Ramadass run');
    ram = run_mode_via_sim_profile('ramadass');
    logmsg('Finished Ramadass run');
    ram_ckpt_file = fullfile(results_dir, ['checkpoint_after_ram_' run_ts '.mat']);
    save(ram_ckpt_file, 'ram', '-v7.3');
    logmsg(['Saved checkpoint: ' ram_ckpt_file]);

    logmsg('Starting Yang run');
    yang = run_mode_via_sim_profile('yang');
    logmsg('Finished Yang run');
    yang_ckpt_file = fullfile(results_dir, ['checkpoint_after_yang_' run_ts '.mat']);
    save(yang_ckpt_file, 'yang', '-v7.3');
    logmsg(['Saved checkpoint: ' yang_ckpt_file]);

    logmsg('Building comparison metrics');

    [t_cmp, ram_j_s, yang_j_s] = overlap_series(ram.time, ram.j_s, yang.time, yang.j_s);
    [~, ram_j_sei, yang_j_sei] = overlap_series(ram.time, ram.j_sei, yang.time, yang.j_sei);
    [~, ram_d_dot, yang_d_dot] = overlap_series(ram.time, ram.delta_sei_dot, yang.time, yang.delta_sei_dot);
    [~, ram_d_dot_num, yang_d_dot_num] = overlap_series(ram.time, ram.delta_sei_dot_num, yang.time, yang.delta_sei_dot_num);
    [~, ram_d, yang_d] = overlap_series(ram.time, ram.delta_sei, yang.time, yang.delta_sei);
    [~, ram_obj_i0s, yang_obj_i0s] = overlap_series(ram.time, ram.obj_p_i0s, yang.time, yang.obj_p_i0s);
    [~, ram_eta_n_side, yang_eta_n_side] = overlap_series(ram.time, ram.eta_n_side, yang.time, yang.eta_n_side);
    [~, ram_phi_se, yang_phi_se] = overlap_series(ram.time, ram.phi_se, yang.time, yang.phi_se);
    [~, ram_eta_s, yang_eta_s] = overlap_series(ram.time, ram.eta_s, yang.time, yang.eta_s);

    logmsg('Generating comparison figures');
    figure('Name', 'SPMe degradation: Ramadass vs Yang', 'Color', 'w');

    subplot(3,2,1);
    plot(t_cmp, ram_j_s, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, yang_j_s, 'r--', 'LineWidth', 1.2);
    plot(t_cmp, yang_j_sei, 'k:', 'LineWidth', 1.2); grid on;
    ylabel('Flux [mol/(m^2 s)]');
    title('Side-reaction fluxes (Ramadass j_s, Yang j_s and j_{sei})');
    legend('Ramadass j_s', 'Yang j_s', 'Yang j_{sei}', 'Location', 'best');

    subplot(3,2,2);
    plot(t_cmp, ram_d_dot, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, yang_d_dot, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('d\delta_{sei}/dt [m/s]');
    title('\delta_{sei} growth rate');
    legend('Ramadass', 'Yang', 'Location', 'best');

    subplot(3,2,3);
    plot(t_cmp, ram_obj_i0s, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, yang_obj_i0s, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('i_{0s} [A/m^2]');
    title('Side-reaction exchange current density parameter obj.p.i0s');
    legend('Ramadass obj.p.i0s', 'Yang obj.p.i0s', 'Location', 'best');

    subplot(3,2,4);
    plot(yang.time, yang.c_s_ec, 'g-', 'LineWidth', 1.2); hold on;
    plot(yang.time, yang.c0_ec*ones(size(yang.time)), 'k--', 'LineWidth', 1.1); grid on;
    ylabel('Concentration');
    title('Yang only: c_{s,ec} vs c_{0,ec}');
    legend('c_{s,ec}', 'c_{0,ec}', 'Location', 'best');

    subplot(3,2,5);
    plot(t_cmp, ram_d*1e9, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, yang_d*1e9, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('\delta_{sei} [nm]');
    title('SEI thickness');
    legend('Ramadass', 'Yang', 'Location', 'best');

    subplot(3,2,6);
    plot(t_cmp, ram_eta_n_side, 'b-', 'LineWidth', 1.2); hold on;
    plot(t_cmp, yang_eta_n_side, 'b--', 'LineWidth', 1.1);
    plot(t_cmp, ram_phi_se, 'r-', 'LineWidth', 1.2);
    plot(t_cmp, yang_phi_se, 'r--', 'LineWidth', 1.1);
    plot(t_cmp, ram_eta_s, 'k-', 'LineWidth', 1.2);
    plot(t_cmp, yang_eta_s, 'k--', 'LineWidth', 1.1); grid on;
    xlabel('Time [s]');
    ylabel('Potential [V]');
    title('\eta_n^{side}, \phi_{se}, \eta_s (Ramadass vs Yang)');
    legend('Ram \eta_n^{side}', 'Yang \eta_n^{side}', ...
           'Ram \phi_{se}', 'Yang \phi_{se}', ...
           'Ram \eta_s', 'Yang \eta_s', 'Location', 'best');

    figure('Name', 'Discharge capacity by cycle (integral method)', 'Color', 'w');
    plot(ram.discharge_cycle_id, ram.discharge_Ah_integral, 'o-b', ...
        'LineWidth', 1.3, 'MarkerSize', 6); hold on;
    plot(yang.discharge_cycle_id, yang.discharge_Ah_integral, 's-r', ...
        'LineWidth', 1.3, 'MarkerSize', 6);
    grid on;
    xlabel('Cycle index');
    ylabel('Discharged capacity [Ah]');
    title('Discharge-only capacity per cycle (integral method)');
    legend('Ramadass', 'Yang', 'Location', 'best');

    out = struct();
    out.ramadass = ram;
    out.yang = yang;
    out.compare.time = t_cmp;
    out.compare.j_s_rmse = sqrt(mean((ram_j_s - yang_j_s).^2));
    out.compare.j_sei_rmse = sqrt(mean((ram_j_sei - yang_j_sei).^2));
    out.compare.delta_sei_dot_rmse = sqrt(mean((ram_d_dot - yang_d_dot).^2));
    out.compare.delta_sei_rmse = sqrt(mean((ram_d - yang_d).^2));
    out.compare.ramadass_dot_vs_num_rmse = sqrt(mean((ram_d_dot - ram_d_dot_num).^2));
    out.compare.yang_dot_vs_num_rmse = sqrt(mean((yang_d_dot - yang_d_dot_num).^2));
    out.compare.discharge_capacity_cycle_rmse = overlap_capacity_rmse( ...
        ram.discharge_cycle_id, ram.discharge_Ah_integral, ...
        yang.discharge_cycle_id, yang.discharge_Ah_integral);

    fprintf('\nDelta-SEI consistency check (analytic vs numerical d(delta)/dt):\n');
    fprintf('  Ramadass RMSE: %.3e [m/s]\n', out.compare.ramadass_dot_vs_num_rmse);
    fprintf('  Yang     RMSE: %.3e [m/s]\n\n', out.compare.yang_dot_vs_num_rmse);
    logmsg('Comparison metrics complete');

    function r = run_mode_via_sim_profile(mode_name)
        logmsg(['Launching sim_profile_sequence_LGM50 for mode=' mode_name]);
        sim_seq = sim_profile_sequence_LGM50('SPMe', init_soc_charge, sequence, ...
            'Tamb', t_amb, 'DeltaT', delta_t, ...
            'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
            'DegradationModel', mode_name, ...
            'ElectrolyteModel', 'fdm', ...
            'SPMeOutputs', {'j_side', 'j_sei', 'delta_sei_dot', 'eta_n_side', 'phi_se', 'eta_s'}, ...
            'Repeat', 1, ...
            'SaveOutput', false);

        seg = post_compute_side_terms(sim_seq, mode_name);
        r = struct();
        r.mode = mode_name;
        r.time = seg.time;
        r.current = seg.current;
        r.voltage = seg.voltage;
        r.j_s = seg.j_s;
        r.j_sei = seg.j_sei;
        r.delta_sei_dot = seg.delta_sei_dot;
        r.delta_sei_dot_num = seg.delta_sei_dot_num;
        r.obj_p_i0s = seg.obj_p_i0s;
        r.eta_n_side = seg.eta_n_side;
        r.phi_se = seg.phi_se;
        r.eta_s = seg.eta_s;
        r.c_s_ec = seg.c_s_ec;
        r.c0_ec = seg.c0_ec;
        r.delta_sei = seg.delta_sei;
        [r.discharge_cycle_id, r.discharge_Ah_integral] = extract_discharge_capacity_integral(sim_seq);
        logmsg(['Completed post-processing for mode=' mode_name]);
    end

    function logmsg(message)
        fprintf('[%s] %s\n', datestr(now, 31), message);
    end
end

function seg = post_compute_side_terms(sim, mode_name)
    nT = numel(sim.time);
    if ~isfield(sim, 'spme_outputs')
        error('sim.spme_outputs is missing. Request SPMeOutputs in sim_profile_sequence_LGM50 call.');
    end
    j_s = sim.spme_outputs.j_side(:);
    j_sei = sim.spme_outputs.j_sei(:);
    dsei_dot = sim.spme_outputs.delta_sei_dot(:);
    dsei_dot_num = zeros(nT,1);
    dsei = sim.x(end, :).';
    obj_p_i0s = zeros(nT,1);
    eta_n_side = sim.spme_outputs.eta_n_side(:);
    phi_se = sim.spme_outputs.phi_se(:);
    eta_s = sim.spme_outputs.eta_s(:);
    c_s_ec = nan(nT,1);
    for k = 1:nT
        if strcmpi(mode_name, 'yang')
            c_s_ec_k = j_sei(k) / sim.p.Faraday * dsei(k) / sim.p.D_ec + sim.p.c0_ec;
            obj_p_i0s(k) = sim.p.k0_sei * c_s_ec_k * sim.p.Faraday;
            c_s_ec(k) = c_s_ec_k;
        else
            obj_p_i0s(k) = sim.p.i0s;
        end
    end
    if nT >= 2
        dt = sim.time(2) - sim.time(1);
        dsei_dot_num(1:end-1) = diff(dsei) / dt;
        dsei_dot_num(end) = dsei_dot_num(end-1);
    end

    seg = struct();
    seg.time = sim.time(:);
    seg.current = sim.current(:);
    seg.voltage = sim.voltage(:);
    seg.j_s = j_s;
    seg.j_sei = j_sei;
    seg.delta_sei_dot = dsei_dot;
    seg.delta_sei_dot_num = dsei_dot_num;
    seg.obj_p_i0s = obj_p_i0s;
    seg.eta_n_side = eta_n_side;
    seg.phi_se = phi_se;
    seg.eta_s = eta_s;
    seg.c_s_ec = c_s_ec;
    seg.c0_ec = sim.p.c0_ec;
    seg.delta_sei = dsei;
end

function [t_cmp, y1_cmp, y2_cmp] = overlap_series(t1, y1, t2, y2)
    t_end = min(t1(end), t2(end));
    t_cmp = t1(t1 <= t_end);
    y1_cmp = y1(1:numel(t_cmp));
    y2_cmp = interp1(t2, y2, t_cmp, 'linear');
end

function [cycle_id, discharge_Ah] = extract_discharge_capacity_integral(sim)
    if ~isfield(sim, 'segment_log') || isempty(sim.segment_log)
        cycle_id = [];
        discharge_Ah = [];
        return;
    end
    seg_log = sim.segment_log;
    n = numel(seg_log);
    cycle_raw = [];
    discharge_raw = [];
    for i = 1:n
        if strcmp(seg_log(i).name, 'CC') && seg_log(i).crate > 0
            if isfield(seg_log, 'capacity_cycle')
                cyc = seg_log(i).capacity_cycle;
            else
                % Fallback for older logs that only have repeat-cycle index.
                cyc = seg_log(i).cycle;
            end
            cycle_raw(end+1,1) = cyc; %#ok<AGROW>
            discharge_raw(end+1,1) = seg_log(i).discharged_Ah; %#ok<AGROW>
        end
    end

    if isempty(cycle_raw)
        cycle_id = [];
        discharge_Ah = [];
        return;
    end

    [cycle_id, ~, ic] = unique(cycle_raw, 'stable');
    discharge_Ah = accumarray(ic, discharge_raw, [], @sum);
end

function e = overlap_capacity_rmse(cyc1, q1, cyc2, q2)
    if isempty(cyc1) || isempty(cyc2)
        e = NaN;
        return;
    end

    [cyc1_u, q1_u] = aggregate_by_cycle(cyc1(:), q1(:));
    [cyc2_u, q2_u] = aggregate_by_cycle(cyc2(:), q2(:));

    cyc_common = intersect(cyc1_u, cyc2_u, 'stable');
    if isempty(cyc_common)
        e = NaN;
        return;
    end
    [~, i1] = ismember(cyc_common, cyc1_u);
    [~, i2] = ismember(cyc_common, cyc2_u);
    q1c = q1_u(i1);
    q2c = q2_u(i2);
    e = sqrt(mean((q1c - q2c).^2, 'omitnan'));
end

function [cyc_u, q_u] = aggregate_by_cycle(cyc, q)
    valid = isfinite(cyc) & isfinite(q);
    cyc = cyc(valid);
    q = q(valid);
    if isempty(cyc)
        cyc_u = [];
        q_u = [];
        return;
    end
    [cyc_u, ~, ic] = unique(cyc, 'stable');
    q_u = accumarray(ic, q, [], @mean);
end

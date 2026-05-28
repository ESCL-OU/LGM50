function out = test_spme_compare(model_a, model_b)
%Compare two SPMe degradation models.

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
    allowed_models = {'ramadass', 'yang', 'safari'};

    [model_a, model_b] = resolve_model_pair(model_a, model_b, allowed_models);

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

    logmsg(['Starting ' model_a ' run']);
    t_model_a = tic;
    res_a = run_mode_via_sim_profile(model_a);
    elapsed_a = toc(t_model_a);
    logmsg(sprintf('Finished %s run (%.3f s)', model_a, elapsed_a));
    ckpt_a_file = fullfile(results_dir, ['checkpoint_after_' model_a '_' run_ts '.mat']);
    save(ckpt_a_file, 'res_a', '-v7.3');
    logmsg(['Saved checkpoint: ' ckpt_a_file]);

    logmsg(['Starting ' model_b ' run']);
    t_model_b = tic;
    res_b = run_mode_via_sim_profile(model_b);
    elapsed_b = toc(t_model_b);
    logmsg(sprintf('Finished %s run (%.3f s)', model_b, elapsed_b));
    ckpt_b_file = fullfile(results_dir, ['checkpoint_after_' model_b '_' run_ts '.mat']);
    save(ckpt_b_file, 'res_b', '-v7.3');
    logmsg(['Saved checkpoint: ' ckpt_b_file]);

    logmsg('Building comparison metrics');

    [t_cmp, a_j_s, b_j_s] = overlap_series(res_a.time, res_a.j_s, res_b.time, res_b.j_s);
    [~, a_j_sei, b_j_sei] = overlap_series(res_a.time, res_a.j_sei, res_b.time, res_b.j_sei);
    [~, a_d_dot, b_d_dot] = overlap_series(res_a.time, res_a.delta_sei_dot, res_b.time, res_b.delta_sei_dot);
    [~, a_d_dot_num, b_d_dot_num] = overlap_series(res_a.time, res_a.delta_sei_dot_num, res_b.time, res_b.delta_sei_dot_num);
    [~, a_d, b_d] = overlap_series(res_a.time, res_a.delta_sei, res_b.time, res_b.delta_sei);
    [~, a_obj_i0s, b_obj_i0s] = overlap_series(res_a.time, res_a.obj_p_i0s, res_b.time, res_b.obj_p_i0s);
    [~, a_eta_n_side, b_eta_n_side] = overlap_series(res_a.time, res_a.eta_n_side, res_b.time, res_b.eta_n_side);
    [~, a_phi_se, b_phi_se] = overlap_series(res_a.time, res_a.phi_se, res_b.time, res_b.phi_se);
    [~, a_eta_s, b_eta_s] = overlap_series(res_a.time, res_a.eta_s, res_b.time, res_b.eta_s);

    logmsg('Generating comparison figures');
    label_a = model_label(model_a);
    label_b = model_label(model_b);
    figure('Name', ['SPMe degradation: ' label_a ' vs ' label_b], 'Color', 'w');

    subplot(3,2,1);
    plot(t_cmp, a_j_s, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, b_j_s, 'r--', 'LineWidth', 1.2);
    plot(t_cmp, b_j_sei, 'k:', 'LineWidth', 1.2); grid on;
    ylabel('Flux [mol/(m^2 s)]');
    title(['Side-reaction fluxes (' label_a ' j_s, ' label_b ' j_s and j_{sei})']);
    legend([label_a ' j_s'], [label_b ' j_s'], [label_b ' j_{sei}'], 'Location', 'best');

    subplot(3,2,2);
    plot(t_cmp, a_d_dot, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, b_d_dot, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('d\delta_{sei}/dt [m/s]');
    title('\delta_{sei} growth rate');
    legend(label_a, label_b, 'Location', 'best');

    subplot(3,2,3);
    plot(t_cmp, a_obj_i0s, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, b_obj_i0s, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('i_{0s} [A/m^2]');
    title('Side-reaction exchange current density parameter obj.p.i0s');
    legend([label_a ' obj.p.i0s'], [label_b ' obj.p.i0s'], 'Location', 'best');

    subplot(3,2,4);
    plot(res_a.time, res_a.c_s_ec, 'g-', 'LineWidth', 1.2); hold on;
    plot(res_b.time, res_b.c_s_ec, 'm--', 'LineWidth', 1.2);
    plot(res_a.time, res_a.c0_ec*ones(size(res_a.time)), 'k:', 'LineWidth', 1.1); grid on;
    ylabel('Concentration');
    title('c_{s,ec} and c_{0,ec}');
    legend([label_a ' c_{s,ec}'], [label_b ' c_{s,ec}'], 'c_{0,ec}', 'Location', 'best');

    subplot(3,2,5);
    plot(t_cmp, a_d*1e9, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_cmp, b_d*1e9, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('\delta_{sei} [nm]');
    title('SEI thickness');
    legend(label_a, label_b, 'Location', 'best');

    subplot(3,2,6);
    plot(t_cmp, a_eta_n_side, 'b-', 'LineWidth', 1.2); hold on;
    plot(t_cmp, b_eta_n_side, 'b--', 'LineWidth', 1.1);
    plot(t_cmp, a_phi_se, 'r-', 'LineWidth', 1.2);
    plot(t_cmp, b_phi_se, 'r--', 'LineWidth', 1.1);
    plot(t_cmp, a_eta_s, 'k-', 'LineWidth', 1.2);
    plot(t_cmp, b_eta_s, 'k--', 'LineWidth', 1.1); grid on;
    xlabel('Time [s]');
    ylabel('Potential [V]');
    title(['\eta_n^{side}, \phi_{se}, \eta_s (' label_a ' vs ' label_b ')']);
    legend([label_a ' \eta_n^{side}'], [label_b ' \eta_n^{side}'], ...
           [label_a ' \phi_{se}'], [label_b ' \phi_{se}'], ...
           [label_a ' \eta_s'], [label_b ' \eta_s'], 'Location', 'best');

    figure('Name', 'Discharge capacity by cycle (integral method)', 'Color', 'w');
    plot(res_a.discharge_cycle_id, res_a.discharge_Ah_integral, 'o-b', ...
        'LineWidth', 1.3, 'MarkerSize', 6); hold on;
    plot(res_b.discharge_cycle_id, res_b.discharge_Ah_integral, 's-r', ...
        'LineWidth', 1.3, 'MarkerSize', 6);
    grid on;
    xlabel('Cycle index');
    ylabel('Discharged capacity [Ah]');
    title('Discharge-only capacity per cycle (integral method)');
    legend(label_a, label_b, 'Location', 'best');

    out = struct();
    out.model_a = model_a;
    out.model_b = model_b;
    out.results.(model_a) = res_a;
    out.results.(model_b) = res_b;
    out.timing.seconds.(model_a) = elapsed_a;
    out.timing.seconds.(model_b) = elapsed_b;
    out.compare.time = t_cmp;
    out.compare.j_s_rmse = sqrt(mean((a_j_s - b_j_s).^2));
    out.compare.j_sei_rmse = sqrt(mean((a_j_sei - b_j_sei).^2));
    out.compare.delta_sei_dot_rmse = sqrt(mean((a_d_dot - b_d_dot).^2));
    out.compare.delta_sei_rmse = sqrt(mean((a_d - b_d).^2));
    out.compare.model_a_dot_vs_num_rmse = sqrt(mean((a_d_dot - a_d_dot_num).^2));
    out.compare.model_b_dot_vs_num_rmse = sqrt(mean((b_d_dot - b_d_dot_num).^2));
    out.compare.discharge_capacity_cycle_rmse = overlap_capacity_rmse( ...
        res_a.discharge_cycle_id, res_a.discharge_Ah_integral, ...
        res_b.discharge_cycle_id, res_b.discharge_Ah_integral);

    fprintf('\nDelta-SEI consistency check (analytic vs numerical d(delta)/dt):\n');
    fprintf('  %-8s RMSE: %.3e [m/s]\n', label_a, out.compare.model_a_dot_vs_num_rmse);
    fprintf('  %-8s RMSE: %.3e [m/s]\n\n', label_b, out.compare.model_b_dot_vs_num_rmse);
    fprintf('Simulation time comparison:\n');
    fprintf('  %-8s : %.3f s\n', label_a, elapsed_a);
    fprintf('  %-8s : %.3f s\n', label_b, elapsed_b);
    fprintf('  %-8s/%-8s speed ratio: %.3f\n\n', label_a, label_b, elapsed_a / elapsed_b);
    logmsg('Comparison metrics complete');

    function r = run_mode_via_sim_profile(mode_name)
        logmsg(['Launching sim_profile_sequence_LGM50 for mode=' mode_name]);
        sim_seq = sim_profile_sequence_LGM50('SPMe', init_soc_charge, sequence, ...
            'Tamb', t_amb, 'DeltaT', delta_t, ...
            'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
            'DegradationModel', mode_name, ...
            'ElectrolyteModel', 'fdm', ...
            'SPMeOutputs', {'j_side', 'j_sei', 'delta_sei', 'delta_sei_dot', 'eta_n_side', 'phi_se', 'eta_s', 'c_s_ec'}, ...
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

    function [m1, m2] = resolve_model_pair(m1, m2, allowed)
        if nargin < 2 || isempty(m1)
            m1 = strtrim(lower(input('Enter first model (ramadass/yang/safari): ', 's')));
        else
            m1 = strtrim(lower(string(m1)));
            m1 = char(m1);
        end
        if nargin < 3 || isempty(m2)
            m2 = strtrim(lower(input('Enter second model (ramadass/yang/safari): ', 's')));
        else
            m2 = strtrim(lower(string(m2)));
            m2 = char(m2);
        end
        validate_model_name(m1, allowed, 'first');
        validate_model_name(m2, allowed, 'second');
        if strcmp(m1, m2)
            error('Please choose two different models. Both were "%s".', m1);
        end
    end

    function validate_model_name(model_name, allowed, position_name)
        if ~any(strcmp(model_name, allowed))
            error('Invalid %s model "%s". Allowed models: %s.', ...
                position_name, model_name, strjoin(allowed, ', '));
        end
    end

    function label = model_label(model_name)
        switch lower(model_name)
            case 'ramadass'
                label = 'Ramadass';
            case 'yang'
                label = 'Yang';
            case 'safari'
                label = 'Safari';
            otherwise
                label = model_name;
        end
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
    if isfield(sim.spme_outputs, 'delta_sei')
        dsei = sim.spme_outputs.delta_sei(:);
    elseif isfield(sim, 'delta_sei')
        dsei = sim.delta_sei(:);
    else
        dsei = sim.x(end, :).';
    end
    obj_p_i0s = zeros(nT,1);
    eta_n_side = sim.spme_outputs.eta_n_side(:);
    phi_se = sim.spme_outputs.phi_se(:);
    eta_s = sim.spme_outputs.eta_s(:);
    if isfield(sim.spme_outputs, 'c_s_ec')
        c_s_ec = sim.spme_outputs.c_s_ec(:);
    else
        c_s_ec = nan(nT,1);
    end
    for k = 1:nT
        if strcmpi(mode_name, 'yang') || strcmpi(mode_name, 'safari')
            c_s_ec_k = c_s_ec(k);
            if ~isfinite(c_s_ec_k)
                c_s_ec_k = j_sei(k) / sim.p.Faraday * dsei(k) / sim.p.D_ec + sim.p.c0_ec;
                c_s_ec(k) = c_s_ec_k;
            end
            if isfield(sim.p, 'k0_sei') && isfinite(c_s_ec_k)
                obj_p_i0s(k) = sim.p.k0_sei * c_s_ec_k * sim.p.Faraday;
            else
                obj_p_i0s(k) = sim.p.i0s;
            end
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

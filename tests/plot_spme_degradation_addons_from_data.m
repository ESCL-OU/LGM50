function out = plot_spme_degradation_addons_from_data(data)
%PLOT_SPME_DEGRADATION_ADDONS_FROM_DATA Visualize existing degradation run data.
%
% This function does not run a simulation. It accepts one of:
%   - output struct from test_spme_degradation_addons_sequence_plot
%   - checkpoint struct saved by sim_profile_sequence_LGM50
%   - sim struct from sim_profile_sequence_LGM50
%   - path to a checkpoint .mat file
%
% Usage after a successful test run:
%   plot_spme_degradation_addons_from_data(ans)
%
% Usage after a failed run with a saved checkpoint:
%   plot_spme_degradation_addons_from_data("results/SPMe/checkpoint_file.mat")

    if nargin < 1 || isempty(data)
        data = evalin('base', 'ans');
    end

    [sim, source_label] = extract_sim(data);
    assert(isstruct(sim), 'Input data does not contain a sim struct.');
    assert(isfield(sim, 'time') && isfield(sim, 'voltage') && isfield(sim, 'soc'), ...
        'sim must contain time, voltage, and soc.');

    t_hr = sim.time(:) / 3600;
    current = get_sim_field(sim, 'current');
    voltage = sim.voltage(:);
    soc = sim.soc(:);

    delta_sei = get_output(sim, 'delta_sei');
    delta_li = get_output(sim, 'delta_li');
    delta_film = get_output(sim, 'delta_film');
    epsilon_e_n = get_output(sim, 'epsilon_e_n');
    j_sei = get_output(sim, 'j_sei');
    j_plating = get_output(sim, 'j_plating');
    j_side = get_output(sim, 'j_side');
    R_film = get_output(sim, 'R_film');

    if isfield(sim, 'capacity_cycle')
        cycle_edges = find(diff(sim.capacity_cycle(:)) > 0) + 1;
        cycle_switch_t_hr = sim.time(cycle_edges) / 3600;
    else
        cycle_switch_t_hr = [];
    end

    [cycle_id, charged_ah, discharged_ah, net_ah] = extract_capacity_log(sim);
    cycle_marker_interval = 100;
    delta_sei_end = values_at_cycle_end(sim, delta_sei, cycle_id);
    delta_li_end = values_at_cycle_end(sim, delta_li, cycle_id);
    epsilon_e_n_end = values_at_cycle_end(sim, epsilon_e_n, cycle_id);
    [n_li_remaining_end, n_li_lost_sei_end, n_li_lost_plating_end, ...
        n_li_total_accounted_end, n_li_initial] = cycle_end_lithium_inventory( ...
        sim, cycle_id, delta_sei, delta_li);
    n_li_solid_end = cycle_end_solid_lithium_inventory(sim, cycle_id);
    init_cs_sample_interval = 500;
    [init_cs_cycle_id, n_li_init_cs_sampled_end, init_cs_sample_table] = sampled_init_cs_lithium_inventory( ...
        sim, cycle_id, init_cs_sample_interval);

    figure('Name', ['SPMe degradation visualization: ' source_label], 'Color', 'w');
    tl = tiledlayout(4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, current, 'm-', 'LineWidth', 1.0); grid on;
    ylabel('Current [A]');
    title('Current');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, voltage, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('Voltage [V]');
    title('Voltage');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, soc, 'b-', 'LineWidth', 1.0); grid on;
    ylabel('SoC [-]');
    title('SoC');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, delta_sei * 1e9, 'r-', 'LineWidth', 1.1); hold on;
    plot(t_hr, delta_li * 1e9, 'b-', 'LineWidth', 1.1);
    plot(t_hr, delta_film * 1e9, 'k--', 'LineWidth', 1.0);
    grid on;
    ylabel('Thickness [nm]');
    legend('\delta_{sei}', '\delta_{Li}', '\delta_{film}', 'Location', 'best');
    title('Film Growth');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, epsilon_e_n, 'g-', 'LineWidth', 1.1); grid on;
    ylabel('\epsilon_{e,n} [-]');
    title('Anode Porosity');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, j_sei, 'r-', 'LineWidth', 1.0); hold on;
    plot(t_hr, j_plating, 'b-', 'LineWidth', 1.0);
    plot(t_hr, j_side, 'k--', 'LineWidth', 1.0);
    grid on;
    ylabel('Flux [mol/(m^2 s)]');
    legend('j_{sei}', 'j_{plating}', 'j_{side}', 'Location', 'best');
    title('Side-Reaction Fluxes');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, R_film, 'Color', [0.5 0 0.5], 'LineWidth', 1.1); grid on;
    xlabel('Time [h]');
    ylabel('R_{film} [\Omega m^2]');
    title('Film Resistance');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(cycle_id, charged_ah, 'o-b', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, discharged_ah, 's-r', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, net_ah, 'd-k', 'LineWidth', 1.0, 'MarkerSize', 4);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Capacity [Ah]');
    legend('Charged', 'Discharged', 'Net', 'Location', 'best');
    title('Capacity');

    title(tl, ['Existing SPMe degradation data: ' source_label], 'Interpreter', 'none');

    figure('Name', ['Cycle-end summary: ' source_label], 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(cycle_id, delta_sei_end * 1e9, 'o-r', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, delta_li_end * 1e9, 's-b', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('End thickness [nm]');
    legend('\delta_{sei}', '\delta_{Li}', 'Location', 'best');

    nexttile;
    plot(cycle_id, epsilon_e_n_end, 'o-g', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('End \epsilon_{e,n} [-]');

    nexttile;
    plot(cycle_id, discharged_ah, 's-r', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, charged_ah, 'o-b', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Capacity [Ah]');
    legend('Discharged', 'Charged', 'Location', 'best');

    nexttile;
    if isempty(discharged_ah)
        capacity_retention = discharged_ah;
    else
        capacity_retention = discharged_ah / max(discharged_ah(1), eps);
    end
    plot(cycle_id, capacity_retention, 'd-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Discharge capacity retention [-]');

    figure('Name', ['Cycle-end lithium inventory: ' source_label], 'Color', 'w');
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(cycle_id, n_li_remaining_end, 'o-g', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, n_li_lost_sei_end, 's-r', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, n_li_lost_plating_end, '^-b', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, n_li_total_accounted_end, 'd-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, n_li_solid_end, 'x-c', 'LineWidth', 1.1, 'MarkerSize', 6);
    plot(init_cs_cycle_id, n_li_init_cs_sampled_end, 'p', ...
        'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
        'LineStyle', 'none', 'MarkerSize', 8);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Lithium inventory [mol]');
    title('Cycle-end inventory');
    legend('Remaining cyclable Li', 'Lost to SEI', 'Lost to Li plating', ...
        'Remaining + losses', 'Solid Li from c_s path', ...
        sprintf('degraded init\\_cs every %d cycles', init_cs_sample_interval), 'Location', 'best');

    nexttile;
    if isfinite(n_li_initial) && n_li_initial > 0
        plot(cycle_id, 100 * n_li_remaining_end / n_li_initial, ...
            'o-g', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
        plot(cycle_id, 100 * n_li_lost_sei_end / n_li_initial, ...
            's-r', 'LineWidth', 1.1, 'MarkerSize', 5);
        plot(cycle_id, 100 * n_li_lost_plating_end / n_li_initial, ...
            '^-b', 'LineWidth', 1.1, 'MarkerSize', 5);
        plot(cycle_id, 100 * n_li_total_accounted_end / n_li_initial, ...
            'd-k', 'LineWidth', 1.1, 'MarkerSize', 5);
        plot(cycle_id, 100 * n_li_solid_end / n_li_initial, ...
            'x-c', 'LineWidth', 1.1, 'MarkerSize', 6);
        plot(init_cs_cycle_id, 100 * n_li_init_cs_sampled_end / n_li_initial, 'p', ...
            'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
            'LineStyle', 'none', 'MarkerSize', 8);
    else
        plot(cycle_id, nan(size(cycle_id)), 'o-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    end
    grid on;
    xlabel('Capacity cycle');
    ylabel('Fraction of initial cyclable Li [%]');
    title(sprintf('Initial cyclable Li = %.6g mol', n_li_initial));
    legend('Remaining cyclable Li', 'Lost to SEI', 'Lost to Li plating', ...
        'Remaining + losses', 'Solid Li from c_s path', ...
        sprintf('degraded init\\_cs every %d cycles', init_cs_sample_interval), 'Location', 'best');

    out = struct();
    out.source_label = source_label;
    out.sim = sim;
    out.cycle_id = cycle_id;
    out.charged_Ah = charged_ah;
    out.discharged_Ah = discharged_ah;
    out.net_Ah = net_ah;
    out.capacity_retention = capacity_retention;
    out.cycle_marker_interval = cycle_marker_interval;
    out.delta_sei = delta_sei;
    out.delta_li = delta_li;
    out.delta_film = delta_film;
    out.epsilon_e_n = epsilon_e_n;
    out.j_sei = j_sei;
    out.j_plating = j_plating;
    out.j_side = j_side;
    out.R_film = R_film;
    out.delta_sei_end = delta_sei_end;
    out.delta_li_end = delta_li_end;
    out.epsilon_e_n_end = epsilon_e_n_end;
    out.n_li_initial = n_li_initial;
    out.n_li_remaining_end = n_li_remaining_end;
    out.n_li_lost_sei_end = n_li_lost_sei_end;
    out.n_li_lost_plating_end = n_li_lost_plating_end;
    out.n_li_total_accounted_end = n_li_total_accounted_end;
    out.n_li_solid_end = n_li_solid_end;
    out.init_cs_sample_interval = init_cs_sample_interval;
    out.init_cs_cycle_id = init_cs_cycle_id;
    out.n_li_init_cs_sampled_end = n_li_init_cs_sampled_end;
    out.init_cs_sample_table = init_cs_sample_table;
end

function [sim, source_label] = extract_sim(data)
    source_label = 'workspace data';
    if ischar(data) || isstring(data)
        source_label = char(data);
        loaded = load(source_label);
        if isfield(loaded, 'checkpoint')
            data = loaded.checkpoint;
        elseif isfield(loaded, 'sim')
            data = loaded.sim;
        else
            names = fieldnames(loaded);
            error('No checkpoint or sim variable found in %s. Variables: %s', ...
                source_label, strjoin(names, ', '));
        end
    end

    if isstruct(data) && isfield(data, 'sim')
        sim = data.sim;
    elseif isstruct(data) && isfield(data, 'checkpoint') && isfield(data.checkpoint, 'sim')
        sim = data.checkpoint.sim;
    elseif isstruct(data) && isfield(data, 'time') && isfield(data, 'voltage')
        sim = data;
    else
        error('Could not find sim data in input.');
    end
end

function y = get_sim_field(sim, name)
    if isfield(sim, name)
        y = sim.(name)(:);
    else
        y = nan(size(sim.time(:)));
    end
    y = align_to_time(y, sim.time);
end

function y = get_output(sim, name)
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, name)
        y = sim.spme_outputs.(name)(:);
    elseif isfield(sim, name)
        y = sim.(name)(:);
    else
        y = nan(size(sim.time(:)));
    end
    y = align_to_time(y, sim.time);
end

function [cycle_id, charged_ah, discharged_ah, net_ah] = extract_capacity_log(sim)
    if isfield(sim, 'capacity_cycle_log') && ~isempty(sim.capacity_cycle_log)
        cycle_log = sim.capacity_cycle_log;
    elseif isfield(sim, 'segment_log') && ~isempty(sim.segment_log)
        cycle_log = build_capacity_cycle_log(sim.segment_log);
    else
        cycle_log = struct('cycle', {}, 'charged_Ah', {}, 'discharged_Ah', {}, 'net_Ah', {});
    end

    n_cycles = numel(cycle_log);
    cycle_id = zeros(n_cycles, 1);
    charged_ah = zeros(n_cycles, 1);
    discharged_ah = zeros(n_cycles, 1);
    net_ah = zeros(n_cycles, 1);
    for k = 1:n_cycles
        cycle_id(k) = cycle_log(k).cycle;
        charged_ah(k) = cycle_log(k).charged_Ah;
        discharged_ah(k) = cycle_log(k).discharged_Ah;
        net_ah(k) = cycle_log(k).net_Ah;
    end
end

function cycle_log = build_capacity_cycle_log(segment_log)
    cycle_ids = [segment_log.capacity_cycle];
    unique_cycles = unique(cycle_ids, 'stable');
    cycle_log = struct('cycle', {}, 'charged_Ah', {}, 'discharged_Ah', {}, 'net_Ah', {});
    for k = 1:numel(unique_cycles)
        cyc = unique_cycles(k);
        idx = find(cycle_ids == cyc);
        charged_ah = sum([segment_log(idx).charged_Ah]);
        discharged_ah = sum([segment_log(idx).discharged_Ah]);
        cycle_log(end+1) = struct( ... %#ok<AGROW>
            'cycle', cyc, ...
            'charged_Ah', charged_ah, ...
            'discharged_Ah', discharged_ah, ...
            'net_Ah', discharged_ah - charged_ah);
    end
end

function y_end = values_at_cycle_end(sim, y, cycle_id)
    y_end = nan(numel(cycle_id), 1);
    if ~isfield(sim, 'capacity_cycle')
        return;
    end
    for k = 1:numel(cycle_id)
        idx_end = find(sim.capacity_cycle == cycle_id(k), 1, 'last');
        if ~isempty(idx_end) && idx_end <= numel(y)
            y_end(k) = y(idx_end);
        end
    end
end

function y = align_to_time(y, t)
    y = y(:);
    n = numel(t);
    if numel(y) == n
        return;
    end
    if isempty(y)
        y = nan(n, 1);
    elseif numel(y) > n
        y = y(1:n);
    else
        y = [y; y(end) * ones(n - numel(y), 1)];
    end
end

function [n_remaining, n_lost_sei, n_lost_plating, n_total_accounted, n_initial] = ...
        cycle_end_lithium_inventory(sim, cycle_id, delta_sei, delta_li)
    if ~isfield(sim, 'p') || ~isstruct(sim.p)
        n_initial = nan;
        n_remaining = nan(size(cycle_id));
        n_lost_sei = nan(size(cycle_id));
        n_lost_plating = nan(size(cycle_id));
        n_total_accounted = nan(size(cycle_id));
        return;
    end

    p = sim.p;
    active_area_n = p.a_s_n * p.Area_n * p.L_n;
    n_initial = initial_cyclable_lithium(sim);

    delta_sei_end = values_at_cycle_end(sim, delta_sei, cycle_id);
    delta_li_end = values_at_cycle_end(sim, delta_li, cycle_id);

    [M_sei, rho_sei] = sei_molar_mass_density(p);
    n_lost_sei = 2 * rho_sei / M_sei * max(delta_sei_end, 0) * active_area_n;
    n_lost_plating = p.rho_li / p.M_li * max(delta_li_end, 0) * active_area_n;

    if isfinite(n_initial) && n_initial > 0
        n_remaining = n_initial - n_lost_sei - n_lost_plating;
        n_total_accounted = n_remaining + n_lost_sei + n_lost_plating;
    else
        n_remaining = nan(size(cycle_id));
        n_total_accounted = nan(size(cycle_id));
    end
end

function n_initial = initial_cyclable_lithium(sim)
    p = sim.p;
    if isfield(p, 'n_Li_s') && isfinite(p.n_Li_s) && p.n_Li_s > 0
        n_initial = p.n_Li_s;
        return;
    end

    if isfield(sim, 'init_csn0') && isfield(sim, 'init_csp0')
        n_initial = sim.init_csn0 * p.epsilon_s_n * p.Area_n * p.L_n + ...
            sim.init_csp0 * p.epsilon_s_p * p.Area_p * p.L_p;
    else
        n_initial = nan;
    end
end

function n_solid_end = cycle_end_solid_lithium_inventory(sim, cycle_id)
    n_solid_end = nan(numel(cycle_id), 1);
    if ~isfield(sim, 'p') || ~isstruct(sim.p)
        return;
    end

    p = sim.p;
    [cs_ave_n, cs_ave_p] = solid_average_concentrations(sim);
    if isempty(cs_ave_n) || isempty(cs_ave_p)
        return;
    end

    n_solid = cs_ave_n(:) * p.epsilon_s_n * p.Area_n * p.L_n + ...
        cs_ave_p(:) * p.epsilon_s_p * p.Area_p * p.L_p;
    n_solid_end = values_at_cycle_end(sim, n_solid, cycle_id);
end

function [sample_cycle_id, n_sampled, sample_table] = sampled_init_cs_lithium_inventory(sim, cycle_id, interval)
    sample_cycle_id = cycle_id(mod(cycle_id, interval) == 0);
    n_sampled = nan(size(sample_cycle_id));
    sample_table = struct('cycle', {}, 'time_s', {}, 'soc', {}, ...
        'delta_sei', {}, 'delta_li', {}, 'epsilon_e_n_state', {}, ...
        'n_lost_sei', {}, 'n_lost_plating', {}, ...
        'csn0_original', {}, 'csp0_original', {}, 'n_Li_s_original', {}, ...
        'csn0_degraded_param', {}, 'csp0_degraded_param', {}, 'n_Li_s_degraded_param', {}, ...
        'veq_degraded', {}, 'csn0_initcs', {}, 'csp0_initcs', {}, ...
        'n_li_initcs', {});
    if isempty(sample_cycle_id) || ~isfield(sim, 'p') || ~isstruct(sim.p) || ~isfield(sim, 'soc')
        return;
    end

    p = sim.p;
    add_param_path_for_init_cs();
    delta_sei = get_output(sim, 'delta_sei');
    delta_li = get_output(sim, 'delta_li');
    epsilon_e_n = get_output(sim, 'epsilon_e_n');
    for k = 1:numel(sample_cycle_id)
        idx_end = find(sim.capacity_cycle == sample_cycle_id(k), 1, 'last');
        if isempty(idx_end) || idx_end > numel(sim.soc)
            continue;
        end
        soc_k = min(max(sim.soc(idx_end), 0), 1);
        try
            dsei_k = value_or_zero(delta_sei, idx_end);
            dli_k = value_or_zero(delta_li, idx_end);
            een_k = value_or_default(epsilon_e_n, idx_end, p.epsilon_e_n);
            [p_deg, n_lost_sei, n_lost_plating] = degraded_params_for_init_cs(p, dsei_k, dli_k);
            [csn0, csp0, veq] = init_cs_from_soc_with_params(p_deg, soc_k);
            n_sampled(k) = csn0 * p_deg.epsilon_s_n * p_deg.Area_n * p_deg.L_n + ...
                csp0 * p_deg.epsilon_s_p * p_deg.Area_p * p_deg.L_p;
            sample_table(k) = struct( ... %#ok<AGROW>
                'cycle', sample_cycle_id(k), ...
                'time_s', sim.time(idx_end), ...
                'soc', soc_k, ...
                'delta_sei', dsei_k, ...
                'delta_li', dli_k, ...
                'epsilon_e_n_state', een_k, ...
                'n_lost_sei', n_lost_sei, ...
                'n_lost_plating', n_lost_plating, ...
                'csn0_original', p.csn0, ...
                'csp0_original', p.csp0, ...
                'n_Li_s_original', get_field_or_nan(p, 'n_Li_s'), ...
                'csn0_degraded_param', p_deg.csn0, ...
                'csp0_degraded_param', p_deg.csp0, ...
                'n_Li_s_degraded_param', get_field_or_nan(p_deg, 'n_Li_s'), ...
                'veq_degraded', veq, ...
                'csn0_initcs', csn0, ...
                'csp0_initcs', csp0, ...
                'n_li_initcs', n_sampled(k));
        catch
            n_sampled(k) = nan;
        end
    end
end

function [p_deg, n_lost_sei, n_lost_plating] = degraded_params_for_init_cs(p, delta_sei, delta_li)
    p_deg = p;
    active_area_n = p.a_s_n * p.Area_n * p.L_n;
    [M_sei, rho_sei] = sei_molar_mass_density(p);
    n_lost_sei = 2 * rho_sei / M_sei * max(delta_sei, 0) * active_area_n;
    n_lost_plating = p.rho_li / p.M_li * max(delta_li, 0) * active_area_n;
    n_lost_total = n_lost_sei + n_lost_plating;

    p_deg.csn0 = p.csn0 - n_lost_total / (p.epsilon_s_n * p.Area_n * p.L_n);
    if isfield(p, 'n_Li_s') && isfinite(p.n_Li_s)
        p_deg.n_Li_s = p.n_Li_s - n_lost_total;
    end
end

function [csn0, csp0, veq] = init_cs_from_soc_with_params(p, soc)
    [cn_low, cp_low] = init_cs(p, p.volt_min);
    [cn_high, cp_high] = init_cs(p, p.volt_max);
    cn_interp = cn_low + soc * (cn_high - cn_low);
    cp_interp = cp_low + soc * (cp_high - cp_low);

    theta_n = cn_interp / p.c_s_n_max;
    theta_p = cp_interp / p.c_s_p_max;
    veq = refPotentialCathode(p, theta_p) - refPotentialAnode(p, theta_n);
    [csn0, csp0] = init_cs(p, veq);
end

function val = value_or_zero(x, idx)
    val = value_or_default(x, idx, 0);
end

function val = value_or_default(x, idx, default_val)
    if ~isempty(x) && idx <= numel(x) && isfinite(x(idx))
        val = x(idx);
    else
        val = default_val;
    end
end

function val = get_field_or_nan(s, name)
    if isstruct(s) && isfield(s, name) && isscalar(s.(name))
        val = s.(name);
    else
        val = nan;
    end
end

function add_param_path_for_init_cs()
    this_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(this_dir);
    param_dir = fullfile(repo_dir, 'param');
    if isfolder(param_dir)
        addpath(param_dir);
    end
end

function [cs_ave_n, cs_ave_p] = solid_average_concentrations(sim)
    cs_ave_n = get_optional_output_raw(sim, 'cs_ave_n');
    cs_ave_p = get_optional_output_raw(sim, 'cs_ave_p');
    if ~isempty(cs_ave_n) && ~isempty(cs_ave_p)
        cs_ave_n = align_to_time(cs_ave_n, sim.time);
        cs_ave_p = align_to_time(cs_ave_p, sim.time);
        return;
    end

    if isfield(sim, 'c_s_n') && isfield(sim, 'c_s_p') && ...
            size(sim.c_s_n, 2) == numel(sim.time) && size(sim.c_s_p, 2) == numel(sim.time)
        cs_ave_n = mean(sim.c_s_n, 1).';
        cs_ave_p = mean(sim.c_s_p, 1).';
        cs_ave_n = align_to_time(cs_ave_n, sim.time);
        cs_ave_p = align_to_time(cs_ave_p, sim.time);
    elseif isfield(sim, 'init_csn0') && isfield(sim, 'init_csp0')
        cs_ave_n = sim.init_csn0 * ones(numel(sim.time), 1);
        cs_ave_p = sim.init_csp0 * ones(numel(sim.time), 1);
    else
        cs_ave_n = [];
        cs_ave_p = [];
    end
end

function y = get_optional_output_raw(sim, name)
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, name)
        y = sim.spme_outputs.(name)(:);
    else
        y = [];
    end
end

function [M_sei, rho_sei] = sei_molar_mass_density(p)
    if isfield(p, 'degradation_model') && strcmpi(char(p.degradation_model), 'ramadass') && ...
            isfield(p, 'M_P') && isfield(p, 'rho_P') && p.M_P > 0 && p.rho_P > 0
        M_sei = p.M_P;
        rho_sei = p.rho_P;
    else
        M_sei = p.M_sei;
        rho_sei = p.rho_sei;
    end
end

function mark_capacity_cycle_edges(cycle_switch_t_hr, marker_interval)
    if nargin < 2 || isempty(marker_interval)
        marker_interval = 100;
    end
    cycle_numbers = (2:(numel(cycle_switch_t_hr) + 1)).';
    keep = mod(cycle_numbers, marker_interval) == 0;
    if ~isempty(keep)
        keep(end) = true;
    end
    for idx = find(keep).'
        xline(cycle_switch_t_hr(idx), '--', sprintf('cycle %d', cycle_numbers(idx)), ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
end

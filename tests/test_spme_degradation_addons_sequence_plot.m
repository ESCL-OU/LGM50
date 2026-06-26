function out = test_spme_degradation_addons_sequence_plot(resume_checkpoint_file, enable_diagnostics, electrolyte_model)
%TEST_SPME_DEGRADATION_ADDONS_SEQUENCE_PLOT
% Simulate a repeated SPMe sequence and visualize SEI, Li plating,
% anode porosity, and capacity evolution.
%
% Usage:
%   out = test_spme_degradation_addons_sequence_plot();
%   out = test_spme_degradation_addons_sequence_plot('/path/to/checkpoint.mat');
%   out = test_spme_degradation_addons_sequence_plot('/path/to/checkpoint.mat', true);
%   out = test_spme_degradation_addons_sequence_plot([], false, 'fvm');

    %% User settings
    degradation_model = 'safari';        % 'ramadass', 'yang', or 'safari'
    enable_lithium_plating = true;
    enable_anode_porosity_loss = true;

    init_soc = 0.1;
    sequence1 = {'REST', 0.25, 'CC', -0.2, 'REST', 0.25, 'CC', 0.2};
    sequence2 = {'CC', -1.0, 'CC', 1.0};
    sequence1_repeat = 1;
    sequence2_repeat = 99;
    repeat_count = 20;

    delta_t = 1;
    t_amb = 25;
    Nr = 15;
    Nxn = 20;
    Nxs = 5;
    Nxp = 20;
    memory_light = true;
    show_progress = true;
    save_checkpoint = true;
    cycle_marker_interval = 100;
    if nargin < 1
        resume_checkpoint_file = '';
    end
    if nargin < 2 || isempty(enable_diagnostics)
        enable_diagnostics = false;
    end
    enable_diagnostics = logical(enable_diagnostics);
    if nargin < 3 || isempty(electrolyte_model)
        electrolyte_model = 'fdm';
    end
    electrolyte_model = lower(char(electrolyte_model));
    if ~ismember(electrolyte_model, {'fdm', 'fvm'})
        error('electrolyte_model must be ''fdm'' or ''fvm''.');
    end

    %% Paths
    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));
    results_dir = fullfile(base_dir, 'results', 'SPMe');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end
    run_ts = datestr(now, 'yyyymmdd_HHMMSS');
    checkpoint_file = '';
    diagnostics_file = '';
    if save_checkpoint
        checkpoint_name = sprintf('spme_degradation_addons_checkpoint_%s_%s.mat', ...
            lower(char(degradation_model)), run_ts);
        checkpoint_file = fullfile(results_dir, checkpoint_name);
        diagnostics_file = fullfile(results_dir, sprintf( ...
            'spme_degradation_addons_diagnostics_%s_%s.mat', ...
            lower(char(degradation_model)), run_ts));
        fprintf('Saving warm-start checkpoints using base path:\n  %s\n', checkpoint_file);
        fprintf('  Each completed repeat writes a separate *_repeat####_segment#####.mat file.\n');
        fprintf('  A manifest is saved as *_manifest.mat.\n');
    end

    %% Run sequence
    block_sequence = [repmat(sequence1, 1, sequence1_repeat), ...
        repmat(sequence2, 1, sequence2_repeat)];

    req_outputs = {'delta_sei', 'delta_sei_dot', ...
        'delta_li', 'delta_li_dot', 'delta_film', 'delta_film_dot', ...
        'epsilon_e_n', 'epsilon_e_n_dot', ...
        'j_sei', 'j_plating', 'j_side', 'R_film'};

    initial_state = struct();
    diagnostic_context = struct('test_name', mfilename, ...
        'run_timestamp', run_ts, ...
        'resume_checkpoint_file', resume_checkpoint_file);
    sequence_to_run = block_sequence;
    repeat_to_run = repeat_count;
    init_soc_run = init_soc;
    resumed_from_checkpoint = false;
    if ~isempty(resume_checkpoint_file)
        ckpt_data = load(resume_checkpoint_file, 'checkpoint');
        if ~isfield(ckpt_data, 'checkpoint')
            error('Checkpoint file does not contain variable "checkpoint": %s', resume_checkpoint_file);
        end
        checkpoint = ckpt_data.checkpoint;
        if ~isfield(checkpoint, 'terminal_state') || ~isfield(checkpoint.terminal_state, 'x')
            error('Checkpoint is missing terminal_state.x: %s', resume_checkpoint_file);
        end
        if ~isfield(checkpoint, 'segments') || ~isfield(checkpoint, 'next_segment_idx')
            error('Checkpoint is missing segment continuation metadata: %s', resume_checkpoint_file);
        end
        if checkpoint.next_segment_idx > numel(checkpoint.segments)
            error('Checkpoint already completed all saved segments: %s', resume_checkpoint_file);
        end
        initial_state = checkpoint.terminal_state;
        init_soc_run = min(max(checkpoint.terminal_soc, 0), 1);
        sequence_to_run = segments_to_profile_sequence(checkpoint.segments(checkpoint.next_segment_idx:end));
        repeat_to_run = 1;
        resumed_from_checkpoint = true;
        diagnostic_context.resume_next_segment_idx = checkpoint.next_segment_idx;
        diagnostic_context.resume_n_segments_total = numel(checkpoint.segments);
        diagnostic_context.resume_terminal_soc = checkpoint.terminal_soc;
        diagnostic_context.resume_terminal_voltage = checkpoint.terminal_voltage;
        fprintf('Resuming from checkpoint:\n  %s\n', resume_checkpoint_file);
        fprintf('  next_segment_idx = %d / %d\n', ...
            checkpoint.next_segment_idx, numel(checkpoint.segments));
        fprintf('  terminal_soc = %.12g, init_soc used for setup = %.12g\n', ...
            checkpoint.terminal_soc, init_soc_run);
    end

    sim = sim_profile_sequence_LGM50('SPMe', init_soc_run, sequence_to_run, ...
        'Repeat', repeat_to_run, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', degradation_model, ...
        'EnableLithiumPlating', enable_lithium_plating, ...
        'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
        'ElectrolyteModel', electrolyte_model, ...
        'SPMeOutputs', req_outputs, ...
        'SaveOutput', false, ...
        'ShowProgress', show_progress, ...
        'MemoryLight', memory_light, ...
        'InitialState', initial_state, ...
        'CheckpointFile', checkpoint_file, ...
        'EnableDiagnostics', enable_diagnostics, ...
        'DiagnosticsFile', diagnostics_file, ...
        'DiagnosticContext', diagnostic_context);

    assert(isfield(sim, 'capacity_cycle_log'), 'Missing sim.capacity_cycle_log.');
    assert(isfield(sim, 'spme_outputs'), 'Missing requested SPMe optional outputs.');

    %% Extract time histories
    t_hr = sim.time(:) / 3600;
    current = sim.current(:);
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

    cycle_log = sim.capacity_cycle_log;
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

    cycle_edges = find(diff(sim.capacity_cycle) > 0) + 1;
    cycle_switch_t_hr = sim.time(cycle_edges) / 3600;
    delta_sei_end = values_at_cycle_end(sim, delta_sei);
    delta_li_end = values_at_cycle_end(sim, delta_li);
    epsilon_e_n_end = values_at_cycle_end(sim, epsilon_e_n);
    [n_li_remaining_end, n_li_lost_sei_end, n_li_lost_plating_end, ...
        n_li_total_accounted_end, n_li_initial] = cycle_end_lithium_inventory( ...
        sim, cycle_id, delta_sei, delta_li, degradation_model);

    %% Time-history figure
    figure_name = sprintf('SPMe degradation add-ons: %s, plating=%d, porosity=%d', ...
        degradation_model, enable_lithium_plating, enable_anode_porosity_loss);
    figure('Name', figure_name, 'Color', 'w');
    tl = tiledlayout(4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, current, 'm-', 'LineWidth', 1.0); grid on;
    ylabel('Current [A]');
    title('Sequence input');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, voltage, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('Voltage [V]');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, soc, 'b-', 'LineWidth', 1.0); grid on;
    ylabel('SoC [-]');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, delta_sei * 1e9, 'r-', 'LineWidth', 1.1); hold on;
    plot(t_hr, delta_li * 1e9, 'b-', 'LineWidth', 1.1);
    plot(t_hr, delta_film * 1e9, 'k--', 'LineWidth', 1.0);
    grid on;
    ylabel('Thickness [nm]');
    legend('\delta_{sei}', '\delta_{Li}', '\delta_{film}', 'Location', 'best');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, epsilon_e_n, 'g-', 'LineWidth', 1.1); grid on;
    ylabel('\epsilon_{e,n} [-]');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, j_sei, 'r-', 'LineWidth', 1.0); hold on;
    plot(t_hr, j_plating, 'b-', 'LineWidth', 1.0);
    plot(t_hr, j_side, 'k--', 'LineWidth', 1.0);
    grid on;
    ylabel('Flux [mol/(m^2 s)]');
    legend('j_{sei}', 'j_{plating}', 'j_{side}', 'Location', 'best');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(t_hr, R_film, 'Color', [0.5 0 0.5], 'LineWidth', 1.1); grid on;
    xlabel('Time [h]');
    ylabel('R_{film} [\Omega m^2]');
    mark_capacity_cycle_edges(cycle_switch_t_hr, cycle_marker_interval);

    nexttile;
    plot(cycle_id, charged_ah, 'o-b', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, discharged_ah, 's-r', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, net_ah, 'd-k', 'LineWidth', 1.0, 'MarkerSize', 4);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Capacity [Ah]');
    legend('Charged', 'Discharged', 'Net', 'Location', 'best');

    title(tl, figure_name, 'Interpreter', 'none');

    %% Cycle-end summary figure
    figure('Name', 'Cycle-end degradation and capacity summary', 'Color', 'w');
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
    if n_cycles > 0
        capacity_retention = discharged_ah / max(discharged_ah(1), eps);
    else
        capacity_retention = discharged_ah;
    end
    plot(cycle_id, capacity_retention, 'd-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Discharge capacity retention [-]');

    %% Cycle-end cyclable-lithium inventory accounting
    figure('Name', 'Cycle-end cyclable lithium inventory accounting', 'Color', 'w');
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(cycle_id, n_li_remaining_end, 'o-g', 'LineWidth', 1.1, 'MarkerSize', 5); hold on;
    plot(cycle_id, n_li_lost_sei_end, 's-r', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, n_li_lost_plating_end, '^-b', 'LineWidth', 1.1, 'MarkerSize', 5);
    plot(cycle_id, n_li_total_accounted_end, 'd-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    grid on;
    xlabel('Capacity cycle');
    ylabel('Lithium inventory [mol]');
    title('Cycle-end inventory');
    legend('Remaining cyclable Li', 'Lost to SEI', 'Lost to Li plating', ...
        'Remaining + losses', 'Location', 'best');

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
    else
        plot(cycle_id, nan(size(cycle_id)), 'o-k', 'LineWidth', 1.1, 'MarkerSize', 5);
    end
    grid on;
    xlabel('Capacity cycle');
    ylabel('Fraction of initial cyclable Li [%]');
    title(sprintf('Initial cyclable Li = %.6g mol', n_li_initial));
    legend('Remaining cyclable Li', 'Lost to SEI', 'Lost to Li plating', ...
        'Remaining + losses', 'Location', 'best');

    %% Output struct
    out = struct();
    out.sim = sim;
    out.settings = struct('degradation_model', degradation_model, ...
        'enable_lithium_plating', enable_lithium_plating, ...
        'enable_anode_porosity_loss', enable_anode_porosity_loss, ...
        'init_soc', init_soc, ...
        'sequence1_repeat', sequence1_repeat, ...
        'sequence2_repeat', sequence2_repeat, ...
        'repeat_count', repeat_count, ...
        'cycle_marker_interval', cycle_marker_interval, ...
        'delta_t', delta_t, ...
        't_amb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'electrolyte_model', electrolyte_model);
    out.checkpoint_file = checkpoint_file;
    out.diagnostics_file = diagnostics_file;
    out.enable_diagnostics = enable_diagnostics;
    out.resume_checkpoint_file = resume_checkpoint_file;
    out.resumed_from_checkpoint = resumed_from_checkpoint;
    out.cycle_id = cycle_id;
    out.charged_Ah = charged_ah;
    out.discharged_Ah = discharged_ah;
    out.net_Ah = net_ah;
    out.capacity_retention = capacity_retention;
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
end

function y = get_output(sim, name)
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, name)
        y = sim.spme_outputs.(name)(:);
    else
        y = nan(size(sim.time(:)));
    end
    y = align_to_time(y, sim.time);
end

function y_end = values_at_cycle_end(sim, y)
    cycle_log = sim.capacity_cycle_log;
    y_end = nan(numel(cycle_log), 1);
    for k = 1:numel(cycle_log)
        idx_end = find(sim.capacity_cycle == cycle_log(k).cycle, 1, 'last');
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
        cycle_end_lithium_inventory(sim, cycle_id, delta_sei, delta_li, degradation_model)
    p = sim.p;
    active_area_n = p.a_s_n * p.Area_n * p.L_n;
    n_initial = initial_cyclable_lithium(sim);

    delta_sei_end = values_at_cycle_end(sim, delta_sei);
    delta_li_end = values_at_cycle_end(sim, delta_li);

    [M_sei, rho_sei] = sei_molar_mass_density(p, degradation_model);
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

function [M_sei, rho_sei] = sei_molar_mass_density(p, degradation_model)
    if strcmpi(char(degradation_model), 'ramadass') && ...
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

function profile_sequence = segments_to_profile_sequence(segments)
    profile_sequence = {};
    for i = 1:numel(segments)
        name = upper(char(segments(i).name));
        if strcmp(name, 'CC')
            profile_sequence = [profile_sequence, {'CC', segments(i).crate}]; %#ok<AGROW>
        elseif strcmp(name, 'REST')
            profile_sequence = [profile_sequence, {'REST', segments(i).rest_hours}]; %#ok<AGROW>
        else
            profile_sequence = [profile_sequence, {name}]; %#ok<AGROW>
        end
    end
end

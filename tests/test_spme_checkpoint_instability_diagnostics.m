function out = test_spme_checkpoint_instability_diagnostics(checkpoint_file)
%TEST_SPME_CHECKPOINT_INSTABILITY_DIAGNOSTICS
% Warm-start from a saved degradation checkpoint and run a short diagnostic
% sequence while recording voltage-term decompositions cycle-by-cycle.
%
% Usage:
%   out = test_spme_checkpoint_instability_diagnostics("results/SPMe/checkpoint.mat");

    %% User settings
    diagnostic_sequence = {'CC', -1.0, 'CC', 1.0};
    diagnostic_repeat_count = 100;
    cycle_duration_long_threshold = 1.10; % flag segments >10% longer than median
    max_overlay_cycles = 100;
    save_diagnostic_output = false;
    use_electrolyte_activity = false;

    if nargin < 1 || isempty(checkpoint_file)
        error('Provide a saved checkpoint .mat file path.');
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

    %% Load checkpoint and settings
    loaded = load(checkpoint_file, 'checkpoint');
    if ~isfield(loaded, 'checkpoint')
        error('Checkpoint file does not contain variable "checkpoint": %s', checkpoint_file);
    end
    checkpoint = loaded.checkpoint;
    if ~isfield(checkpoint, 'terminal_state') || ~isfield(checkpoint.terminal_state, 'x')
        error('Checkpoint is missing terminal_state.x.');
    end
    opt = checkpoint.options;

    model_name = 'SPMe';
    if isfield(checkpoint, 'model_name') && ~isempty(checkpoint.model_name)
        model_name = char(checkpoint.model_name);
    end
    degradation_model = char_or_default(get_option(opt, 'DegradationModel', 'safari'), 'safari');
    electrolyte_model = char_or_default(get_option(opt, 'ElectrolyteModel', 'fdm'), 'fdm');
    enable_lithium_plating = get_option(opt, 'EnableLithiumPlating', true);
    enable_anode_porosity_loss = get_option(opt, 'EnableAnodePorosityLoss', true);
    delta_t = get_option(opt, 'DeltaT', 1);
    t_amb = get_option(opt, 'Tamb', 25);
    Nr = get_option(opt, 'Nr', 15);
    Nxn = get_option(opt, 'Nxn', 30);
    Nxs = get_option(opt, 'Nxs', 5);
    Nxp = get_option(opt, 'Nxp', 30);

    req_outputs = {'delta_sei', 'delta_sei_dot', ...
        'delta_li', 'delta_li_dot', 'delta_film', 'delta_film_dot', ...
        'epsilon_e_n', 'epsilon_e_n_dot', ...
        'j_sei', 'j_plating', 'j_side', 'R_film', ...
        'eta_n_side', 'phi_se', 'eta_s'};

    segments = repeat_segments(parse_profile_sequence(diagnostic_sequence), diagnostic_repeat_count);
    n_segments = numel(segments);

    current_state = checkpoint.terminal_state;
    if ~isfield(current_state, 'p') || isempty(current_state.p)
        if isfield(checkpoint, 'last_segment') && isfield(checkpoint.last_segment, 'p')
            current_state.p = checkpoint.last_segment.p;
        elseif isfield(checkpoint, 'sim') && isfield(checkpoint.sim, 'p')
            current_state.p = checkpoint.sim.p;
        else
            error('Checkpoint must contain terminal_state.p, last_segment.p, or sim.p.');
        end
    end
    raw_init_soc = checkpoint.terminal_soc;
    global_time_offset = 0;
    diag = initialize_diag();
    failure = struct('occurred', false, 'segment_idx', NaN, 'message', '', ...
        'raw_init_soc', NaN, 'segment', struct());

    fprintf('Warm-starting diagnostic from checkpoint:\n  %s\n', checkpoint_file);
    fprintf('  terminal_soc = %.12g, terminal_voltage = %.12g\n', ...
        checkpoint.terminal_soc, checkpoint.terminal_voltage);
    fprintf('  running %d diagnostic segments (%d repeats)\n', n_segments, diagnostic_repeat_count);

    %% Run manually segment-by-segment so partial results remain available
    for s = 1:n_segments
        seg = segments(s);
        raw_soc_invalid = ~isscalar(raw_init_soc) || ~isfinite(raw_init_soc) || ...
            raw_init_soc < 0 || raw_init_soc > 1;

        try
            seg_sim = sim_profile_LGM50(model_name, raw_init_soc, segment_to_profile_def(seg), ...
                'DeltaT', delta_t, ...
                'Tamb', t_amb, ...
                'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
                'DegradationModel', degradation_model, ...
                'EnableLithiumPlating', enable_lithium_plating, ...
                'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
                'ElectrolyteModel', electrolyte_model, ...
                'UseElectrolyteActivity', use_electrolyte_activity, ...
                'SPMeOutputs', req_outputs, ...
                'InitialState', current_state, ...
                'SaveOutput', false, ...
                'SuppressSummary', true);
        catch ME
            failure.occurred = true;
            failure.segment_idx = s;
            failure.message = ME.message;
            failure.raw_init_soc = raw_init_soc;
            failure.segment = seg;
            warning('Diagnostic stopped at segment %d: %s', s, ME.message);
            break;
        end

        terms = compute_voltage_terms(seg_sim);
        diag = append_segment_diag(diag, seg_sim, terms, seg, s, ...
            global_time_offset, raw_init_soc, raw_init_soc, raw_soc_invalid);

        current_state = struct('x', seg_sim.x(:, end), 'p', seg_sim.p);
        raw_init_soc = seg_sim.soc(end);
        global_time_offset = diag.time(end);
    end

    diag.failure = failure;
    diag.checkpoint_file = checkpoint_file;
    diag.settings = struct('diagnostic_sequence', {diagnostic_sequence}, ...
        'diagnostic_repeat_count', diagnostic_repeat_count, ...
        'cycle_duration_long_threshold', cycle_duration_long_threshold, ...
        'degradation_model', degradation_model, ...
        'electrolyte_model', electrolyte_model, ...
        'use_electrolyte_activity', use_electrolyte_activity, ...
        'enable_lithium_plating', enable_lithium_plating, ...
        'enable_anode_porosity_loss', enable_anode_porosity_loss, ...
        'delta_t', delta_t, 't_amb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp);

    diag = summarize_durations(diag, cycle_duration_long_threshold);
    plot_duration_summary(diag);
    plot_voltage_term_overlay(diag, 'charge', max_overlay_cycles);
    plot_voltage_term_overlay(diag, 'discharge', max_overlay_cycles);
    plot_degradation_summary(diag);

    out = diag;
    if save_diagnostic_output
        run_ts = datestr(now, 'yyyymmdd_HHMMSS');
        out_file = fullfile(results_dir, ['spme_checkpoint_instability_diagnostics_' run_ts '.mat']);
        save(out_file, 'out', '-v7.3');
        out.saved_file = out_file;
        fprintf('Saved diagnostic output:\n  %s\n', out_file);
    end
end

function diag = initialize_diag()
    diag = struct();
    diag.time = [];
    diag.segment_elapsed = [];
    diag.current = [];
    diag.voltage = [];
    diag.soc = [];
    diag.segment_idx = [];
    diag.diagnostic_cycle = [];
    diag.segment_mode = {};
    diag.raw_init_soc = [];
    diag.setup_soc = [];
    diag.raw_soc_invalid = [];
    diag.ocv = [];
    diag.reaction = [];
    diag.solid_film_ohmic = [];
    diag.electrolyte_cond = [];
    diag.electrolyte_polar = [];
    diag.reconstructed_voltage = [];
    diag.delta_sei = [];
    diag.delta_li = [];
    diag.delta_film = [];
    diag.epsilon_e_n = [];
    diag.R_film = [];
    diag.segment_log = struct('segment_idx', {}, 'diagnostic_cycle', {}, ...
        'mode', {}, 'duration_s', {}, 'raw_init_soc', {}, 'setup_soc', {}, ...
        'raw_soc_invalid', {}, 'start_voltage', {}, 'end_voltage', {}, ...
        'start_soc', {}, 'end_soc', {});
end

function diag = append_segment_diag(diag, seg_sim, terms, seg, segment_idx, time_offset, raw_init_soc, setup_soc, raw_soc_invalid)
    n = numel(seg_sim.time);
    idx = 1:n;
    if ~isempty(diag.time)
        idx = 2:n;
    end
    t_global = time_offset + seg_sim.time(idx);
    elapsed = seg_sim.time(idx) - seg_sim.time(1);
    mode = segment_mode(seg);
    diagnostic_cycle = seg.cycle;

    diag.time = [diag.time; t_global(:)];
    diag.segment_elapsed = [diag.segment_elapsed; elapsed(:)];
    diag.current = [diag.current; seg_sim.current(idx)];
    diag.voltage = [diag.voltage; seg_sim.voltage(idx)];
    diag.soc = [diag.soc; seg_sim.soc(idx)];
    diag.segment_idx = [diag.segment_idx; segment_idx * ones(numel(idx), 1)];
    diag.diagnostic_cycle = [diag.diagnostic_cycle; diagnostic_cycle * ones(numel(idx), 1)];
    diag.segment_mode = [diag.segment_mode; repmat({mode}, numel(idx), 1)];
    diag.raw_init_soc = [diag.raw_init_soc; raw_init_soc * ones(numel(idx), 1)];
    diag.setup_soc = [diag.setup_soc; setup_soc * ones(numel(idx), 1)];
    diag.raw_soc_invalid = [diag.raw_soc_invalid; raw_soc_invalid * ones(numel(idx), 1)];
    diag.ocv = [diag.ocv; terms.ocv(idx)];
    diag.reaction = [diag.reaction; terms.reaction(idx)];
    diag.solid_film_ohmic = [diag.solid_film_ohmic; terms.solid_film_ohmic(idx)];
    diag.electrolyte_cond = [diag.electrolyte_cond; terms.electrolyte_cond(idx)];
    diag.electrolyte_polar = [diag.electrolyte_polar; terms.electrolyte_polar(idx)];
    diag.reconstructed_voltage = [diag.reconstructed_voltage; terms.reconstructed_voltage(idx)];
    diag.delta_sei = [diag.delta_sei; get_output(seg_sim, 'delta_sei', idx)];
    diag.delta_li = [diag.delta_li; get_output(seg_sim, 'delta_li', idx)];
    diag.delta_film = [diag.delta_film; get_output(seg_sim, 'delta_film', idx)];
    diag.epsilon_e_n = [diag.epsilon_e_n; get_output(seg_sim, 'epsilon_e_n', idx)];
    diag.R_film = [diag.R_film; get_output(seg_sim, 'R_film', idx)];

    duration_s = seg_sim.time(end) - seg_sim.time(1);
    diag.segment_log(end+1) = struct( ... %#ok<AGROW>
        'segment_idx', segment_idx, ...
        'diagnostic_cycle', diagnostic_cycle, ...
        'mode', mode, ...
        'duration_s', duration_s, ...
        'raw_init_soc', raw_init_soc, ...
        'setup_soc', setup_soc, ...
        'raw_soc_invalid', raw_soc_invalid, ...
        'start_voltage', seg_sim.voltage(1), ...
        'end_voltage', seg_sim.voltage(end), ...
        'start_soc', seg_sim.soc(1), ...
        'end_soc', seg_sim.soc(end));
end

function terms = compute_voltage_terms(seg_sim)
    if strcmpi(char(seg_sim.p.degradation_model), 'safari')
        solver = 'cvodes';
    else
        solver = 'idas';
    end
    spme = SPMe(seg_sim.p, solver);
    n = numel(seg_sim.time);
    fields = {'ocv', 'reaction', 'solid_film_ohmic', 'electrolyte_cond', ...
        'electrolyte_polar', 'reconstructed_voltage'};
    for f = 1:numel(fields)
        terms.(fields{f}) = zeros(n, 1);
    end
    for k = 1:n
        [eta_n, eta_p, Unref, Upref, V_noVCE, V_electrolyteCond, V_electrolytePolar] = ...
            spme.potentials(seg_sim.x(:, k), seg_sim.current(k));
        terms.ocv(k) = Upref - Unref;
        terms.reaction(k) = eta_p - eta_n;
        terms.solid_film_ohmic(k) = V_noVCE - terms.ocv(k) - terms.reaction(k);
        terms.electrolyte_cond(k) = V_electrolyteCond;
        terms.electrolyte_polar(k) = V_electrolytePolar;
        terms.reconstructed_voltage(k) = V_noVCE + V_electrolyteCond + V_electrolytePolar;
    end
end

function diag = summarize_durations(diag, threshold)
    modes = {'charge', 'discharge'};
    for m = 1:numel(modes)
        mode = modes{m};
        idx = strcmp({diag.segment_log.mode}, mode);
        durations = [diag.segment_log(idx).duration_s].';
        cycles = [diag.segment_log(idx).diagnostic_cycle].';
        if isempty(durations)
            median_duration = NaN;
            long_cycles = [];
        else
            median_duration = median(durations, 'omitnan');
            long_cycles = cycles(durations > threshold * median_duration);
        end
        diag.duration_summary.(mode) = struct('cycle', cycles, ...
            'duration_s', durations, ...
            'median_duration_s', median_duration, ...
            'long_cycles', long_cycles);
    end
end

function plot_duration_summary(diag)
    figure('Name', 'Checkpoint diagnostic: cycle duration and handoff SOC', 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot_duration_mode(diag, 'charge');
    title('Charge Duration');

    nexttile;
    plot_duration_mode(diag, 'discharge');
    title('Discharge Duration');

    nexttile;
    semilogy(max(diag.delta_film, eps) * 1e9, 'k-', 'LineWidth', 1.0); grid on;
    xlabel('Diagnostic sample');
    ylabel('\delta_{film} [nm]');
    title('Film Thickness');

    nexttile;
    plot([diag.segment_log.diagnostic_cycle], [diag.segment_log.raw_init_soc], 'o-k', ...
        'LineWidth', 1.0, 'MarkerSize', 4); hold on;
    bad = [diag.segment_log.raw_soc_invalid];
    if any(bad)
        plot([diag.segment_log(bad).diagnostic_cycle], [diag.segment_log(bad).raw_init_soc], ...
            'xr', 'LineWidth', 1.5, 'MarkerSize', 8);
    end
    yline(0, '--');
    yline(1, '--');
    grid on;
    xlabel('Diagnostic cycle');
    ylabel('Raw handoff SoC [-]');
    title('Segment Handoff SoC');
end

function plot_duration_mode(diag, mode)
    s = diag.duration_summary.(mode);
    plot(s.cycle, s.duration_s / 60, 'o-', 'LineWidth', 1.0, 'MarkerSize', 4); hold on;
    if isfinite(s.median_duration_s)
        yline(s.median_duration_s / 60, '--k', 'median');
    end
    if ~isempty(s.long_cycles)
        is_long = ismember(s.cycle, s.long_cycles);
        plot(s.cycle(is_long), s.duration_s(is_long) / 60, 'xr', ...
            'LineWidth', 1.5, 'MarkerSize', 8);
    end
    grid on;
    xlabel('Diagnostic cycle');
    ylabel('Duration [min]');
end

function plot_voltage_term_overlay(diag, mode, max_cycles)
    idx_mode = strcmp(diag.segment_mode, mode);
    cycles = unique(diag.diagnostic_cycle(idx_mode), 'stable');
    if numel(cycles) > max_cycles
        cycles = cycles(end-max_cycles+1:end);
    end
    if isempty(cycles)
        return;
    end

    fig = figure('Name', ['Checkpoint diagnostic voltage terms: ' mode], 'Color', 'w');
    tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    colors = parula(numel(cycles));
    plot_term('voltage', 'Terminal V [V]');
    plot_term('ocv', 'OCV [V]');
    plot_term('reaction', '\eta_p-\eta_n [V]');
    plot_term('solid_film_ohmic', 'Solid/film ohmic [V]');
    plot_term('electrolyte_cond', 'Electrolyte cond [V]');
    plot_term('electrolyte_polar', 'Electrolyte polar [V]');
    cb = colorbar;
    cb.Layout.Tile = 'east';
    colormap(fig, colors);
    if numel(cycles) > 1
        caxis([cycles(1), cycles(end)]);
    else
        caxis([cycles(1)-0.5, cycles(1)+0.5]);
    end
    cb.Label.String = 'Diagnostic cycle';
    title(tl, sprintf('%s voltage terms (dark/blue older, yellow newer)', mode), ...
        'Interpreter', 'none');

    function plot_term(field_name, y_label)
        nexttile;
        hold on;
        for ii = 1:numel(cycles)
            cyc = cycles(ii);
            idx = idx_mode & diag.diagnostic_cycle == cyc;
            x = diag.segment_elapsed(idx) / 60;
            y = diag.(field_name)(idx);
            is_long = ismember(cyc, diag.duration_summary.(mode).long_cycles);
            if is_long
                plot(x, y, 'Color', colors(ii, :), 'LineWidth', 1.8);
            else
                plot(x, y, 'Color', colors(ii, :), 'LineWidth', 0.9);
            end
        end
        grid on;
        xlabel('Segment elapsed time [min]');
        ylabel(y_label);
        title(y_label);
    end
end

function plot_degradation_summary(diag)
    figure('Name', 'Checkpoint diagnostic degradation states', 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(diag.time / 3600, diag.delta_sei * 1e9, 'r-', 'LineWidth', 1.0); hold on;
    plot(diag.time / 3600, diag.delta_li * 1e9, 'b-', 'LineWidth', 1.0);
    plot(diag.time / 3600, diag.delta_film * 1e9, 'k--', 'LineWidth', 1.0);
    grid on;
    xlabel('Diagnostic time [h]');
    ylabel('Thickness [nm]');
    legend('\delta_{sei}', '\delta_{Li}', '\delta_{film}', 'Location', 'best');

    nexttile;
    plot(diag.time / 3600, diag.epsilon_e_n, 'g-', 'LineWidth', 1.0);
    grid on;
    xlabel('Diagnostic time [h]');
    ylabel('\epsilon_{e,n} [-]');

    nexttile;
    plot(diag.time / 3600, diag.R_film, 'Color', [0.5 0 0.5], 'LineWidth', 1.0);
    grid on;
    xlabel('Diagnostic time [h]');
    ylabel('R_{film} [\Omega m^2]');

    nexttile;
    plot(diag.time / 3600, diag.voltage - diag.reconstructed_voltage, 'k-', 'LineWidth', 1.0);
    grid on;
    xlabel('Diagnostic time [h]');
    ylabel('V - V_{recon} [V]');
end

function y = get_output(seg_sim, name, idx)
    if isfield(seg_sim, 'spme_outputs') && isfield(seg_sim.spme_outputs, name)
        v = seg_sim.spme_outputs.(name)(:);
    elseif isfield(seg_sim, name)
        v = seg_sim.(name)(:);
    else
        v = nan(numel(seg_sim.time), 1);
    end
    y = v(idx);
end

function segments = repeat_segments(segments_base, repeat_count)
    segments = struct('name', {}, 'crate', {}, 'rest_hours', {}, 'cycle', {});
    for c = 1:repeat_count
        for s = 1:numel(segments_base)
            segments(end+1) = struct('name', segments_base(s).name, ... %#ok<AGROW>
                'crate', segments_base(s).crate, ...
                'rest_hours', segments_base(s).rest_hours, ...
                'cycle', c);
        end
    end
end

function segments = parse_profile_sequence(profile_sequence)
    segments = struct('name', {}, 'crate', {}, 'rest_hours', {});
    i = 1;
    while i <= numel(profile_sequence)
        token = upper(char(profile_sequence{i}));
        if strcmp(token, 'CC')
            segments(end+1) = struct('name', 'CC', ... %#ok<AGROW>
                'crate', profile_sequence{i+1}, ...
                'rest_hours', NaN);
            i = i + 2;
        elseif strcmp(token, 'REST')
            segments(end+1) = struct('name', 'REST', ... %#ok<AGROW>
                'crate', 0, ...
                'rest_hours', profile_sequence{i+1});
            i = i + 2;
        else
            segments(end+1) = struct('name', token, ... %#ok<AGROW>
                'crate', NaN, ...
                'rest_hours', NaN);
            i = i + 1;
        end
    end
end

function prof = segment_to_profile_def(seg)
    if strcmp(seg.name, 'CC')
        prof = {'CC', seg.crate};
    elseif strcmp(seg.name, 'REST')
        prof = {'REST', seg.rest_hours};
    else
        prof = seg.name;
    end
end

function mode = segment_mode(seg)
    if strcmp(seg.name, 'CC') && seg.crate < 0
        mode = 'charge';
    elseif strcmp(seg.name, 'CC') && seg.crate > 0
        mode = 'discharge';
    elseif strcmp(seg.name, 'REST')
        mode = 'rest';
    else
        mode = 'dynamic';
    end
end

function x = get_option(opt, name, default_value)
    if isstruct(opt) && isfield(opt, name) && ~isempty(opt.(name))
        x = opt.(name);
    else
        x = default_value;
    end
end

function y = char_or_default(x, default_value)
    if isempty(x)
        y = default_value;
    else
        y = char(x);
    end
end


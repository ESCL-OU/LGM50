function out = test_spme_warmstart_failure_trace_plot(checkpoint_file, manifest_file)
%TEST_SPME_WARMSTART_FAILURE_TRACE_PLOT Replay a warm start and plot failure states.
%
% Usage:
%   out = test_spme_warmstart_failure_trace_plot("results/SPMe/checkpoint.mat");
%   out = test_spme_warmstart_failure_trace_plot("results/SPMe/checkpoint.mat", ...
%       "results/SPMe/checkpoint_manifest.mat");
%   out = test_spme_warmstart_failure_trace_plot("results/SPMe/checkpoint_manifest.mat");
%
% The replay uses the saved checkpoint continuation state and the same
% remaining current-profile segments. It stops at the solver failure, then
% plots the stitched successful history plus the partial failing segment.

    if nargin < 1 || isempty(checkpoint_file)
        error('Provide a warm-start checkpoint .mat file path.');
    end
    if nargin < 2
        manifest_file = '';
    end

    ensure_test_paths();
    if is_manifest_file(checkpoint_file)
        out = replay_manifest_checkpoints_one_cycle(checkpoint_file);
        return;
    end

    results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'SPMe');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    loaded = load(checkpoint_file, 'checkpoint');
    if ~isfield(loaded, 'checkpoint')
        error('Checkpoint file does not contain variable "checkpoint": %s', checkpoint_file);
    end
    checkpoint = loaded.checkpoint;
    if ~isfield(checkpoint, 'terminal_state') || ~isfield(checkpoint.terminal_state, 'x')
        error('Checkpoint is missing terminal_state.x.');
    end
    if ~isfield(checkpoint, 'segments') || ~isfield(checkpoint, 'next_segment_idx')
        error('Checkpoint is missing segments/next_segment_idx metadata.');
    end

    opt = checkpoint.options;
    model_name = char_or_default(get_field_or(checkpoint, 'model_name', 'SPMe'), 'SPMe');
    degradation_model = char_or_default(get_option(opt, 'DegradationModel', 'safari'), 'safari');
    electrolyte_model = char_or_default(get_option(opt, 'ElectrolyteModel', 'fdm'), 'fdm');
    enable_lithium_plating = get_option(opt, 'EnableLithiumPlating', true);
    enable_anode_porosity_loss = get_option(opt, 'EnableAnodePorosityLoss', true);
    delta_t = get_option(opt, 'DeltaT', 1);
    t_amb = get_option(opt, 'Tamb', 25);
    Nr = get_option(opt, 'Nr', 15);
    Nxn = get_option(opt, 'Nxn', 20);
    Nxs = get_option(opt, 'Nxs', 5);
    Nxp = get_option(opt, 'Nxp', 20);

    req_outputs = default_req_outputs();

    run_ts = datestr(now, 'yyyymmdd_HHMMSS');
    diagnostics_file = fullfile(results_dir, ...
        ['spme_warmstart_failure_trace_' run_ts '_step.mat']);

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

    remaining_segments = checkpoint.segments(checkpoint.next_segment_idx:end);
    init_soc = min(max(checkpoint.terminal_soc, 0), 1);
    trace = initialize_trace();
    failure = struct('occurred', false, 'voltage_terminated', false, ...
        'local_segment_idx', NaN, 'original_segment_idx', NaN, ...
        'message', '', 'diagnostics_file', diagnostics_file);
    global_time_offset = 0;

    fprintf('Replaying warm start from checkpoint:\n  %s\n', checkpoint_file);
    fprintf('  first original segment = %d / %d\n', ...
        checkpoint.next_segment_idx, numel(checkpoint.segments));

    for s = 1:numel(remaining_segments)
        seg = remaining_segments(s);
        original_segment_idx = checkpoint.next_segment_idx + s - 1;
        context = struct('test_name', mfilename, ...
            'checkpoint_file', checkpoint_file, ...
            'local_segment_idx', s, ...
            'original_segment_idx', original_segment_idx, ...
            'segment', seg);

        try
            seg_sim = sim_profile_LGM50(model_name, init_soc, segment_to_profile_def(seg), ...
                'DeltaT', delta_t, ...
                'Tamb', t_amb, ...
                'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
                'DegradationModel', degradation_model, ...
                'EnableLithiumPlating', enable_lithium_plating, ...
                'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
                'ElectrolyteModel', electrolyte_model, ...
                'SPMeOutputs', req_outputs, ...
                'InitialState', current_state, ...
                'DiagnosticsFile', diagnostics_file, ...
                'DiagnosticContext', context, ...
                'SaveOutput', false, ...
                'SuppressSummary', true);
        catch ME
            failure.occurred = true;
            failure.local_segment_idx = s;
            failure.original_segment_idx = original_segment_idx;
            failure.message = ME.message;
            failure.segment = seg;
            if isfile(diagnostics_file)
                loaded_diag = load(diagnostics_file, 'diagnostic');
                if isfield(loaded_diag, 'diagnostic') && isfield(loaded_diag.diagnostic, 'partial_history')
                    trace = append_partial_history(trace, loaded_diag.diagnostic.partial_history, ...
                        seg, s, original_segment_idx, global_time_offset);
                    failure.step_idx = loaded_diag.diagnostic.step_idx;
                    failure.step_time_s = loaded_diag.diagnostic.step_time_s;
                    failure.x0_summary = loaded_diag.diagnostic.x0_summary;
                end
            end
            warning('Replay stopped at original segment %d: %s', original_segment_idx, ME.message);
            break;
        end

        trace = append_segment_sim(trace, seg_sim, seg, s, original_segment_idx, global_time_offset);
        current_state = struct('x', seg_sim.x(:, end), 'p', seg_sim.p);
        init_soc = seg_sim.soc(end);
        global_time_offset = trace.time(end);
    end

    trace.failure = failure;
    trace.checkpoint_file = checkpoint_file;
    trace.diagnostics_file = diagnostics_file;
    plot_trace(trace);
    print_nan_origin_summary(trace);
    checkpoint_overlay = load_checkpoint_manifest_overlay(checkpoint, checkpoint_file, manifest_file);
    if ~isempty(checkpoint_overlay)
        plot_checkpoint_manifest_overlay(checkpoint_overlay, trace);
    end
    trace.checkpoint_overlay = checkpoint_overlay;
    out = trace;
end

function trace = initialize_trace()
    trace = struct();
    trace.time = [];
    trace.segment_elapsed = [];
    trace.current = [];
    trace.voltage = [];
    trace.soc = [];
    trace.local_segment_idx = [];
    trace.original_segment_idx = [];
    trace.segment_mode = {};
    trace.epsilon_e_n_output = [];
    trace.epsilon_e_n_state = [];
    trace.delta_sei = [];
    trace.delta_li = [];
    trace.delta_film = [];
    trace.R_film = [];
    trace.c_e_min = [];
    trace.c_e_mean = [];
    trace.c_e_max = [];
    trace.ce0n = [];
    trace.cens = [];
    trace.cesp = [];
    trace.ce0p = [];
    trace.electrolyte_log_arg_min = [];
    trace.c_ec_sei_min = [];
    trace.c_ec_sei_mean = [];
    trace.c_ec_sei_max = [];
    trace.V_electrolytePolar = [];
    trace.cen_bar = [];
    trace.ces_bar = [];
    trace.cep_bar = [];
end

function trace = append_segment_sim(trace, sim, seg, local_segment_idx, original_segment_idx, time_offset)
    idx = 1:numel(sim.time);
    if ~isempty(trace.time)
        idx = 2:numel(sim.time);
    end
    state = state_series_from_sim(sim);
    trace = append_trace_arrays(trace, sim.time(idx), sim.current(idx), sim.voltage(idx), ...
        sim.soc(idx), state, idx, get_outputs_from_sim(sim), seg, ...
        local_segment_idx, original_segment_idx, time_offset);
end

function trace = append_partial_history(trace, hist, seg, local_segment_idx, original_segment_idx, time_offset)
    idx = 1:numel(hist.time_s);
    if ~isempty(trace.time)
        idx = 2:numel(hist.time_s);
    end
    outputs = get_outputs_from_partial(hist);
    state = hist;
    trace = append_trace_arrays(trace, hist.time_s(idx), hist.current_A(idx), ...
        hist.voltage(idx), hist.soc(idx), state, idx, outputs, seg, ...
        local_segment_idx, original_segment_idx, time_offset);
end

function trace = append_trace_arrays(trace, t_local, current, voltage, soc, state, idx, outputs, ...
        seg, local_segment_idx, original_segment_idx, time_offset)
    n = numel(idx);
    trace.time = [trace.time; time_offset + t_local(:)];
    trace.segment_elapsed = [trace.segment_elapsed; t_local(:)];
    trace.current = [trace.current; current(:)];
    trace.voltage = [trace.voltage; voltage(:)];
    trace.soc = [trace.soc; soc(:)];
    trace.local_segment_idx = [trace.local_segment_idx; local_segment_idx * ones(n, 1)];
    trace.original_segment_idx = [trace.original_segment_idx; original_segment_idx * ones(n, 1)];
    trace.segment_mode = [trace.segment_mode; repmat({segment_mode(seg)}, n, 1)];
    trace.epsilon_e_n_output = [trace.epsilon_e_n_output; get_output_array(outputs, 'epsilon_e_n', idx, n)];
    trace.epsilon_e_n_state = [trace.epsilon_e_n_state; get_state_array(state, 'epsilon_e_n_state', idx, n)];
    trace.delta_sei = [trace.delta_sei; prefer_output_or_state(outputs, state, 'delta_sei', 'delta_sei_state', idx, n)];
    trace.delta_li = [trace.delta_li; prefer_output_or_state(outputs, state, 'delta_li', 'delta_li_state', idx, n)];
    trace.delta_film = [trace.delta_film; get_output_array(outputs, 'delta_film', idx, n)];
    trace.R_film = [trace.R_film; get_output_array(outputs, 'R_film', idx, n)];
    trace.c_e_min = [trace.c_e_min; get_state_array(state, 'c_e_min', idx, n)];
    trace.c_e_mean = [trace.c_e_mean; get_state_array(state, 'c_e_mean', idx, n)];
    trace.c_e_max = [trace.c_e_max; get_state_array(state, 'c_e_max', idx, n)];
    ce0n = prefer_output_or_state(outputs, state, 'ce0n', 'ce0n', idx, n);
    cens = prefer_output_or_state(outputs, state, 'cens', 'cens', idx, n);
    cesp = prefer_output_or_state(outputs, state, 'cesp', 'cesp', idx, n);
    ce0p = prefer_output_or_state(outputs, state, 'ce0p', 'ce0p', idx, n);
    trace.ce0n = [trace.ce0n; ce0n];
    trace.cens = [trace.cens; cens];
    trace.cesp = [trace.cesp; cesp];
    trace.ce0p = [trace.ce0p; ce0p];
    trace.electrolyte_log_arg_min = [trace.electrolyte_log_arg_min; ...
        min([ce0n, cens, cesp, ce0p], [], 2)];
    trace.c_ec_sei_min = [trace.c_ec_sei_min; get_state_array(state, 'c_ec_sei_min', idx, n)];
    trace.c_ec_sei_mean = [trace.c_ec_sei_mean; get_state_array(state, 'c_ec_sei_mean', idx, n)];
    trace.c_ec_sei_max = [trace.c_ec_sei_max; get_state_array(state, 'c_ec_sei_max', idx, n)];
    trace.V_electrolytePolar = [trace.V_electrolytePolar; get_output_array(outputs, 'V_electrolytePolar', idx, n)];
    trace.cen_bar = [trace.cen_bar; get_output_array(outputs, 'cen_bar', idx, n)];
    trace.ces_bar = [trace.ces_bar; get_output_array(outputs, 'ces_bar', idx, n)];
    trace.cep_bar = [trace.cep_bar; get_output_array(outputs, 'cep_bar', idx, n)];
end

function state = state_series_from_sim(sim)
    p = sim.p;
    state = struct();
    c_e_hist = electrolyte_history_from_sim(sim);
    state.c_e_min = min(c_e_hist, [], 1).';
    state.c_e_mean = mean(c_e_hist, 1).';
    state.c_e_max = max(c_e_hist, [], 1).';
    boundary = electrolyte_boundary_history(c_e_hist, p);
    state.ce0n = boundary.ce0n;
    state.cens = boundary.cens;
    state.cesp = boundary.cesp;
    state.ce0p = boundary.ce0p;
    state.electrolyte_log_arg_min = min([state.ce0n, state.cens, state.cesp, state.ce0p], [], 2);
    if isfield(sim, 'epsilon_e_n')
        state.epsilon_e_n_state = sim.epsilon_e_n(:);
    else
        state.epsilon_e_n_state = p.epsilon_e_n * ones(numel(sim.time), 1);
    end
    if isfield(sim, 'delta_sei')
        state.delta_sei_state = sim.delta_sei(:);
    end
    if isfield(sim, 'delta_li')
        state.delta_li_state = sim.delta_li(:);
    end
    if isfield(sim, 'c_ec_sei')
        state.c_ec_sei_min = min(sim.c_ec_sei, [], 1).';
        state.c_ec_sei_mean = mean(sim.c_ec_sei, 1).';
        state.c_ec_sei_max = max(sim.c_ec_sei, [], 1).';
    end
end

function c_e_hist = electrolyte_history_from_sim(sim)
    p = sim.p;
    if isfield(sim, 'x') && ~isempty(sim.x)
        idx_ce_start = 2 * p.Nr - 1;
        idx_ce_end = idx_ce_start + p.Nx - 4;
        if idx_ce_end <= size(sim.x, 1)
            c_e_hist = sim.x(idx_ce_start:idx_ce_end, :);
            return;
        end
    end
    if isfield(sim, 'c_e') && ~isempty(sim.c_e)
        c_e_hist = sim.c_e;
    else
        c_e_hist = nan(p.Nx - 3, numel(sim.time));
    end
end

function boundary = electrolyte_boundary_history(c_e_hist, p)
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, C_ce] = c_e_mats(p);
    n_internal = size(C_ce, 2);
    n_full = n_internal + 4;
    if size(c_e_hist, 1) ~= n_internal && size(c_e_hist, 2) == n_internal
        c_e_hist = c_e_hist.';
    end
    if size(c_e_hist, 1) ~= n_full && size(c_e_hist, 2) == n_full
        c_e_hist = c_e_hist.';
    end
    n = size(c_e_hist, 2);
    boundary = struct();
    boundary.ce0n = nan(n, 1);
    boundary.cens = nan(n, 1);
    boundary.cesp = nan(n, 1);
    boundary.ce0p = nan(n, 1);
    for j = 1:n
        c_e = c_e_hist(:, j);
        if numel(c_e) == n_internal
            c_e_bcs = full(C_ce * c_e);
            c_en = c_e(1:(p.Nxn-1));
            c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
            c_ep = c_e((p.Nxn-1)+p.Nxs:end);
            c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];
        elseif numel(c_e) == n_full
            c_ex = c_e;
        else
            continue;
        end
        boundary.ce0n(j) = c_ex(1);
        boundary.cens(j) = c_ex(p.Nxn+1);
        boundary.cesp(j) = c_ex(p.Nxn+p.Nxs+1);
        boundary.ce0p(j) = c_ex(end);
    end
end

function outputs = get_outputs_from_sim(sim)
    if isfield(sim, 'spme_outputs')
        outputs = sim.spme_outputs;
    else
        outputs = struct();
    end
end

function outputs = get_outputs_from_partial(hist)
    if isfield(hist, 'spme_outputs')
        outputs = hist.spme_outputs;
    else
        outputs = struct();
    end
end

function y = get_output_array(outputs, name, idx, n)
    if isfield(outputs, name)
        v = outputs.(name)(:);
        y = v(idx);
    else
        y = nan(n, 1);
    end
end

function y = get_state_array(state, name, idx, n)
    if isfield(state, name)
        v = state.(name)(:);
        y = v(idx);
    else
        y = nan(n, 1);
    end
end

function y = prefer_output_or_state(outputs, state, output_name, state_name, idx, n)
    y = get_output_array(outputs, output_name, idx, n);
    missing = isnan(y);
    if any(missing)
        ys = get_state_array(state, state_name, idx, n);
        y(missing) = ys(missing);
    end
end

function overlay = load_checkpoint_manifest_overlay(reference_checkpoint, checkpoint_file, manifest_file)
    overlay = [];
    if nargin < 3 || isempty(manifest_file)
        manifest_file = '';
    else
        manifest_file = char(manifest_file);
    end
    if isempty(manifest_file) && isfield(reference_checkpoint, 'checkpoint_manifest_file') && ...
            ~isempty(reference_checkpoint.checkpoint_manifest_file)
        manifest_file = reference_checkpoint.checkpoint_manifest_file;
    elseif isempty(manifest_file)
        manifest_file = checkpoint_manifest_path_from_checkpoint(checkpoint_file);
    end
    if isempty(manifest_file) || ~isfile(manifest_file)
        fprintf('Checkpoint manifest not found; skipping checkpoint overlay.\n');
        return;
    end
    loaded = load(manifest_file, 'checkpoint_manifest');
    if ~isfield(loaded, 'checkpoint_manifest') || ~isfield(loaded.checkpoint_manifest, 'files')
        fprintf('Checkpoint manifest has no file list; skipping checkpoint overlay.\n');
        return;
    end

    files = loaded.checkpoint_manifest.files(:);
    overlay = struct();
    overlay.manifest_file = manifest_file;
    overlay.files = files;
    overlay.segment_idx = [];
    overlay.repeat_idx = [];
    overlay.terminal_time = [];
    overlay.terminal_voltage = [];
    overlay.terminal_soc = [];
    overlay.traces = struct('file', {}, 'segment_idx', {}, 'repeat_idx', {}, ...
        'elapsed_hr', {}, 'current', {}, 'voltage', {}, 'ce0n', {}, 'cens', {}, ...
        'cesp', {}, 'ce0p', {}, 'log_arg_min', {}, 'c_e_min', {}, ...
        'epsilon_e_n', {});

    for i = 1:numel(files)
        file_i = files{i};
        if ~isfile(file_i)
            continue;
        end
        data = load(file_i, 'checkpoint');
        if ~isfield(data, 'checkpoint') || ~isfield(data.checkpoint, 'last_segment')
            continue;
        end
        ckpt = data.checkpoint;
        seg_sim = ckpt.last_segment;
        if ~isstruct(seg_sim) || ~isfield(seg_sim, 'time') || ~isfield(seg_sim, 'p')
            continue;
        end
        tr = checkpoint_segment_trace(seg_sim);
        tr.file = file_i;
        tr.segment_idx = get_field_or(ckpt, 'segment_idx', NaN);
        if isfield(ckpt, 'segments') && isfinite(tr.segment_idx) && ...
                tr.segment_idx >= 1 && tr.segment_idx <= numel(ckpt.segments)
            tr.repeat_idx = ckpt.segments(tr.segment_idx).cycle;
        else
            tr.repeat_idx = NaN;
        end
        tr = orderfields(tr, overlay.traces);
        overlay.traces(end+1) = tr; %#ok<AGROW>
        overlay.segment_idx(end+1, 1) = tr.segment_idx; %#ok<AGROW>
        overlay.repeat_idx(end+1, 1) = tr.repeat_idx; %#ok<AGROW>
        overlay.terminal_time(end+1, 1) = get_field_or(ckpt, 'terminal_time', NaN); %#ok<AGROW>
        overlay.terminal_voltage(end+1, 1) = get_field_or(ckpt, 'terminal_voltage', NaN); %#ok<AGROW>
        overlay.terminal_soc(end+1, 1) = get_field_or(ckpt, 'terminal_soc', NaN); %#ok<AGROW>
    end

    if isempty(overlay.traces)
        overlay = [];
        fprintf('No checkpoint last_segment trajectories found in manifest.\n');
    end
end

function tr = checkpoint_segment_trace(sim)
    outputs = get_outputs_from_sim(sim);
    state = state_series_from_sim(sim);
    idx = (1:numel(sim.time)).';
    n = numel(idx);
    ce0n = prefer_output_or_state(outputs, state, 'ce0n', 'ce0n', idx, n);
    cens = prefer_output_or_state(outputs, state, 'cens', 'cens', idx, n);
    cesp = prefer_output_or_state(outputs, state, 'cesp', 'cesp', idx, n);
    ce0p = prefer_output_or_state(outputs, state, 'ce0p', 'ce0p', idx, n);
    tr = struct();
    tr.elapsed_hr = (sim.time(:) - sim.time(1)) / 3600;
    tr.current = sim.current(:);
    tr.voltage = sim.voltage(:);
    tr.ce0n = ce0n;
    tr.cens = cens;
    tr.cesp = cesp;
    tr.ce0p = ce0p;
    tr.log_arg_min = min([ce0n, cens, cesp, ce0p], [], 2);
    tr.c_e_min = get_state_array(state, 'c_e_min', idx, n);
    tr.epsilon_e_n = prefer_output_or_state(outputs, state, ...
        'epsilon_e_n', 'epsilon_e_n_state', idx, n);
end

function out = plot_manifest_checkpoint_sims(manifest_file)
    manifest_file = char(manifest_file);
    loaded = load(manifest_file, 'checkpoint_manifest');
    if ~isfield(loaded, 'checkpoint_manifest') || ~isfield(loaded.checkpoint_manifest, 'files')
        error('Manifest file does not contain checkpoint_manifest.files: %s', manifest_file);
    end
    files = loaded.checkpoint_manifest.files(:);
    traces = struct('file', {}, 'segment_idx', {}, 'repeat_idx', {}, ...
        'time_hr', {}, 'current', {}, 'voltage', {}, 'soc', {}, ...
        'delta_sei', {}, 'delta_li', {}, 'delta_film', {}, ...
        'epsilon_e_n', {}, 'R_film', {}, 'ce0n', {}, 'cens', {}, ...
        'cesp', {}, 'ce0p', {}, 'log_arg_min', {}, 'c_e_min', {}, ...
        'c_e_mean', {}, 'c_e_max', {});
    segment_idx = [];
    repeat_idx = [];
    terminal_time = [];
    terminal_soc = [];
    terminal_voltage = [];

    fprintf('Loading checkpoint sims from manifest:\n  %s\n', manifest_file);
    for i = 1:numel(files)
        file_i = files{i};
        if ~isfile(file_i)
            warning('Checkpoint file not found, skipping: %s', file_i);
            continue;
        end
        data = load(file_i, 'checkpoint');
        if ~isfield(data, 'checkpoint') || ~isfield(data.checkpoint, 'sim')
            warning('Checkpoint has no checkpoint.sim, skipping: %s', file_i);
            continue;
        end
        ckpt = data.checkpoint;
        sim = ckpt.sim;
        if ~isstruct(sim) || ~isfield(sim, 'time') || ~isfield(sim, 'voltage')
            warning('checkpoint.sim is missing time/voltage, skipping: %s', file_i);
            continue;
        end
        tr = checkpoint_sim_trace(sim);
        tr.file = file_i;
        tr.segment_idx = get_field_or(ckpt, 'segment_idx', NaN);
        if isfield(ckpt, 'segments') && isfinite(tr.segment_idx) && ...
                tr.segment_idx >= 1 && tr.segment_idx <= numel(ckpt.segments)
            tr.repeat_idx = ckpt.segments(tr.segment_idx).cycle;
        else
            tr.repeat_idx = NaN;
        end
        tr = orderfields(tr, traces);
        traces(end+1) = tr; %#ok<AGROW>
        segment_idx(end+1, 1) = tr.segment_idx; %#ok<AGROW>
        repeat_idx(end+1, 1) = tr.repeat_idx; %#ok<AGROW>
        terminal_time(end+1, 1) = get_field_or(ckpt, 'terminal_time', last_or_nan(sim.time)); %#ok<AGROW>
        terminal_soc(end+1, 1) = get_field_or(ckpt, 'terminal_soc', last_or_nan(sim.soc)); %#ok<AGROW>
        terminal_voltage(end+1, 1) = get_field_or(ckpt, 'terminal_voltage', last_or_nan(sim.voltage)); %#ok<AGROW>
    end

    if isempty(traces)
        error('No usable checkpoint.sim histories found in manifest: %s', manifest_file);
    end

    out = struct();
    out.manifest_file = manifest_file;
    out.files = files;
    out.traces = traces;
    out.segment_idx = segment_idx;
    out.repeat_idx = repeat_idx;
    out.terminal_time = terminal_time;
    out.terminal_soc = terminal_soc;
    out.terminal_voltage = terminal_voltage;

    plot_manifest_full_sim_overlay(out);
    last_segment_overlay = load_checkpoint_manifest_overlay(struct(), '', manifest_file);
    if ~isempty(last_segment_overlay)
        plot_checkpoint_manifest_overlay(last_segment_overlay, struct());
        out.last_segment_overlay = last_segment_overlay;
    end
    fprintf('Loaded %d checkpoint sim histories from manifest.\n', numel(traces));
end

function out = replay_manifest_checkpoints_one_cycle(manifest_file)
    manifest_file = char(manifest_file);
    loaded = load(manifest_file, 'checkpoint_manifest');
    if ~isfield(loaded, 'checkpoint_manifest') || ~isfield(loaded.checkpoint_manifest, 'files')
        error('Manifest file does not contain checkpoint_manifest.files: %s', manifest_file);
    end
    files = loaded.checkpoint_manifest.files(:);
    req_outputs = default_req_outputs();
    traces = struct('file', {}, 'segment_idx', {}, 'repeat_idx', {}, ...
        'time_hr', {}, 'current', {}, 'voltage', {}, 'soc', {}, ...
        'delta_sei', {}, 'delta_li', {}, 'delta_film', {}, ...
        'epsilon_e_n', {}, 'R_film', {}, 'ce0n', {}, 'cens', {}, ...
        'cesp', {}, 'ce0p', {}, 'log_arg_min', {}, 'c_e_min', {}, ...
        'c_e_mean', {}, 'c_e_max', {}, 'failure', {});
    segment_idx = [];
    repeat_idx = [];

    fprintf('Warm-start replaying one cycle from each checkpoint in manifest:\n  %s\n', manifest_file);
    for i = 1:numel(files)
        file_i = files{i};
        if ~isfile(file_i)
            warning('Checkpoint file not found, skipping: %s', file_i);
            continue;
        end
        data = load(file_i, 'checkpoint');
        if ~isfield(data, 'checkpoint')
            warning('File has no checkpoint variable, skipping: %s', file_i);
            continue;
        end
        ckpt = data.checkpoint;
        try
            tr = replay_one_cycle_from_checkpoint(ckpt, file_i, req_outputs);
        catch ME
            warning('One-cycle replay failed before producing trace for %s: %s', file_i, ME.message);
            continue;
        end
        tr = orderfields(tr, traces);
        traces(end+1) = tr; %#ok<AGROW>
        segment_idx(end+1, 1) = tr.segment_idx; %#ok<AGROW>
        repeat_idx(end+1, 1) = tr.repeat_idx; %#ok<AGROW>
        fprintf('  replayed checkpoint %d/%d: segment %g, %d samples\n', ...
            i, numel(files), tr.segment_idx, numel(tr.time_hr));
    end

    if isempty(traces)
        error('No one-cycle warm-start replays succeeded from manifest: %s', manifest_file);
    end

    out = struct();
    out.manifest_file = manifest_file;
    out.files = files;
    out.traces = traces;
    out.segment_idx = segment_idx;
    out.repeat_idx = repeat_idx;
    plot_manifest_one_cycle_replay_overlay(out);
end

function tr = replay_one_cycle_from_checkpoint(checkpoint, checkpoint_file, req_outputs)
    if ~isfield(checkpoint, 'terminal_state') || ~isfield(checkpoint.terminal_state, 'x')
        error('Checkpoint is missing terminal_state.x.');
    end
    if ~isfield(checkpoint, 'segments') || ~isfield(checkpoint, 'next_segment_idx')
        error('Checkpoint is missing segments/next_segment_idx.');
    end
    opt = checkpoint.options;
    model_name = char_or_default(get_field_or(checkpoint, 'model_name', 'SPMe'), 'SPMe');
    degradation_model = char_or_default(get_option(opt, 'DegradationModel', 'safari'), 'safari');
    electrolyte_model = char_or_default(get_option(opt, 'ElectrolyteModel', 'fdm'), 'fdm');
    enable_lithium_plating = get_option(opt, 'EnableLithiumPlating', true);
    enable_anode_porosity_loss = get_option(opt, 'EnableAnodePorosityLoss', true);
    delta_t = get_option(opt, 'DeltaT', 1);
    t_amb = get_option(opt, 'Tamb', 25);
    Nr = get_option(opt, 'Nr', 15);
    Nxn = get_option(opt, 'Nxn', 20);
    Nxs = get_option(opt, 'Nxs', 5);
    Nxp = get_option(opt, 'Nxp', 20);

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

    segments = one_capacity_cycle_segments(checkpoint.segments, checkpoint.next_segment_idx);
    trace = initialize_trace();
    init_soc = min(max(checkpoint.terminal_soc, 0), 1);
    time_offset = 0;
    failure = struct('occurred', false, 'message', '', 'segment', struct());
    first_segment_idx = checkpoint.next_segment_idx;

    for s = 1:numel(segments)
        seg = segments(s);
        original_segment_idx = first_segment_idx + s - 1;
        try
            seg_sim = sim_profile_LGM50(model_name, init_soc, segment_to_profile_def(seg), ...
                'DeltaT', delta_t, ...
                'Tamb', t_amb, ...
                'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
                'DegradationModel', degradation_model, ...
                'EnableLithiumPlating', enable_lithium_plating, ...
                'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
                'ElectrolyteModel', electrolyte_model, ...
                'SPMeOutputs', req_outputs, ...
                'InitialState', current_state, ...
                'SaveOutput', false, ...
                'SuppressSummary', true);
        catch ME
            failure.occurred = true;
            failure.message = ME.message;
            failure.segment = seg;
            break;
        end
        trace = append_segment_sim(trace, seg_sim, seg, s, original_segment_idx, time_offset);
        current_state = struct('x', seg_sim.x(:, end), 'p', seg_sim.p);
        init_soc = seg_sim.soc(end);
        time_offset = trace.time(end);
    end

    tr = trace_to_overlay_trace(trace);
    tr.file = checkpoint_file;
    tr.segment_idx = get_field_or(checkpoint, 'segment_idx', first_segment_idx - 1);
    if isfield(checkpoint, 'segments') && isfinite(tr.segment_idx) && ...
            tr.segment_idx >= 1 && tr.segment_idx <= numel(checkpoint.segments)
        tr.repeat_idx = checkpoint.segments(tr.segment_idx).cycle;
    else
        tr.repeat_idx = NaN;
    end
    tr.failure = failure;
end

function segments_one = one_capacity_cycle_segments(segments, start_idx)
    if start_idx > numel(segments)
        segments_one = segments([]);
        return;
    end
    remaining = segments(start_idx:end);
    cap_cycle = capacity_cycle_index_for_segments(remaining);
    keep = cap_cycle == cap_cycle(1);
    last_keep = find(keep, 1, 'last');
    segments_one = remaining(1:last_keep);
end

function cap_cycle = capacity_cycle_index_for_segments(segments)
    n = numel(segments);
    seg_dir = zeros(n, 1);
    for i = 1:n
        if strcmp(segments(i).name, 'CC')
            seg_dir(i) = sign(segments(i).crate);
        elseif strcmp(segments(i).name, 'REST')
            seg_dir(i) = 0;
        else
            seg_dir(i) = 1;
        end
    end
    half_cycle = ones(n, 1);
    last_nonzero_dir = seg_dir(1);
    for i = 2:n
        half_cycle(i) = half_cycle(i-1);
        if seg_dir(i) ~= 0
            if last_nonzero_dir ~= 0 && seg_dir(i) ~= last_nonzero_dir
                half_cycle(i) = half_cycle(i) + 1;
            end
            last_nonzero_dir = seg_dir(i);
        end
    end
    cap_cycle = ceil(half_cycle / 2);
end

function tr = trace_to_overlay_trace(trace)
    tr = struct();
    tr.file = '';
    tr.segment_idx = NaN;
    tr.repeat_idx = NaN;
    tr.time_hr = trace.time(:) / 3600;
    tr.current = trace.current(:);
    tr.voltage = trace.voltage(:);
    tr.soc = trace.soc(:);
    tr.delta_sei = trace.delta_sei(:);
    tr.delta_li = trace.delta_li(:);
    tr.delta_film = trace.delta_film(:);
    tr.epsilon_e_n = trace.epsilon_e_n_output(:);
    tr.R_film = trace.R_film(:);
    tr.ce0n = trace.ce0n(:);
    tr.cens = trace.cens(:);
    tr.cesp = trace.cesp(:);
    tr.ce0p = trace.ce0p(:);
    tr.log_arg_min = trace.electrolyte_log_arg_min(:);
    tr.c_e_min = trace.c_e_min(:);
    tr.c_e_mean = trace.c_e_mean(:);
    tr.c_e_max = trace.c_e_max(:);
    tr.failure = struct();
end

function plot_manifest_one_cycle_replay_overlay(out)
    plot_manifest_full_sim_overlay(out);
    failed = arrayfun(@(x) isfield(x.failure, 'occurred') && x.failure.occurred, out.traces);
    if any(failed)
        fprintf('One-cycle replay failures occurred for %d/%d checkpoints.\n', ...
            nnz(failed), numel(out.traces));
    end
end

function tr = checkpoint_sim_trace(sim)
    outputs = get_outputs_from_sim(sim);
    state = state_series_from_sim(sim);
    idx = (1:numel(sim.time)).';
    n = numel(idx);
    ce0n = prefer_output_or_state(outputs, state, 'ce0n', 'ce0n', idx, n);
    cens = prefer_output_or_state(outputs, state, 'cens', 'cens', idx, n);
    cesp = prefer_output_or_state(outputs, state, 'cesp', 'cesp', idx, n);
    ce0p = prefer_output_or_state(outputs, state, 'ce0p', 'ce0p', idx, n);
    tr = struct();
    tr.file = '';
    tr.segment_idx = NaN;
    tr.repeat_idx = NaN;
    tr.time_hr = sim.time(:) / 3600;
    tr.current = get_sim_array(sim, 'current', n);
    tr.voltage = sim.voltage(:);
    tr.soc = sim.soc(:);
    tr.delta_sei = prefer_output_or_state(outputs, state, 'delta_sei', 'delta_sei_state', idx, n);
    tr.delta_li = prefer_output_or_state(outputs, state, 'delta_li', 'delta_li_state', idx, n);
    tr.delta_film = get_output_array(outputs, 'delta_film', idx, n);
    tr.epsilon_e_n = prefer_output_or_state(outputs, state, ...
        'epsilon_e_n', 'epsilon_e_n_state', idx, n);
    tr.R_film = get_output_array(outputs, 'R_film', idx, n);
    tr.ce0n = ce0n;
    tr.cens = cens;
    tr.cesp = cesp;
    tr.ce0p = ce0p;
    tr.log_arg_min = min([ce0n, cens, cesp, ce0p], [], 2);
    tr.c_e_min = get_state_array(state, 'c_e_min', idx, n);
    tr.c_e_mean = get_state_array(state, 'c_e_mean', idx, n);
    tr.c_e_max = get_state_array(state, 'c_e_max', idx, n);
end

function y = get_sim_array(sim, name, n)
    if isfield(sim, name)
        y = sim.(name)(:);
    else
        y = nan(n, 1);
    end
end

function plot_manifest_full_sim_overlay(out)
    n = numel(out.traces);
    colors = parula(max(n, 2));

    figure('Name', 'Manifest overlay: degradation histories from checkpoint sims', 'Color', 'w');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    plot_overlay_panel(out.traces, colors, 'current', 'Current [A]', 'Current');
    plot_overlay_panel(out.traces, colors, 'voltage', 'Voltage [V]', 'Voltage');
    plot_overlay_panel(out.traces, colors, 'soc', 'SoC [-]', 'SoC');
    plot_overlay_panel(out.traces, colors, 'epsilon_e_n', '\epsilon_{e,n} [-]', 'Anode Porosity');

    nexttile;
    hold on;
    for i = 1:n
        tr = out.traces(i);
        plot(tr.time_hr, tr.delta_sei * 1e9, '-', 'Color', colors(i, :), 'LineWidth', 0.8);
        plot(tr.time_hr, tr.delta_li * 1e9, '--', 'Color', colors(i, :), 'LineWidth', 0.8);
    end
    grid on;
    xlabel('Checkpoint sim time [h]');
    ylabel('Thickness [nm]');
    title('\delta_{sei} solid, \delta_{Li} dashed');

    plot_overlay_panel(out.traces, colors, 'R_film', 'R_{film} [\Omega m^2]', 'Film Resistance');
    add_checkpoint_colorbar(colors, n);

    figure('Name', 'Manifest overlay: electrolyte depletion from checkpoint sims', 'Color', 'w');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_overlay_panel(out.traces, colors, 'ce0n', 'ce0n [mol/m^3]', 'Anode boundary electrolyte');
    plot_overlay_panel(out.traces, colors, 'log_arg_min', 'min log arg [mol/m^3]', 'Minimum log argument');
    plot_overlay_panel(out.traces, colors, 'c_e_min', 'min c_e [mol/m^3]', 'Internal electrolyte minimum');
    plot_overlay_panel(out.traces, colors, 'c_e_mean', 'mean c_e [mol/m^3]', 'Internal electrolyte mean');
    plot_overlay_panel(out.traces, colors, 'c_e_max', 'max c_e [mol/m^3]', 'Internal electrolyte max');
    plot_overlay_panel(out.traces, colors, 'epsilon_e_n', '\epsilon_{e,n} [-]', 'Anode porosity');
    add_checkpoint_colorbar(colors, n);

    figure('Name', 'Manifest overlay: terminal trends from checkpoint sims', 'Color', 'w');
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    x = out.segment_idx(:);
    plot_terminal_panel(x, out.traces, 'ce0n', 'terminal ce0n');
    plot_terminal_panel(x, out.traces, 'log_arg_min', 'terminal min log arg');
    plot_terminal_panel(x, out.traces, 'c_e_min', 'terminal min c_e');
    plot_terminal_panel(x, out.traces, 'c_e_mean', 'terminal mean c_e');
    plot_terminal_panel(x, out.traces, 'epsilon_e_n', 'terminal \epsilon_{e,n}');
    plot(x, out.terminal_voltage, 'o-k', 'LineWidth', 1.0); grid on;
    xlabel('Checkpoint segment index');
    ylabel('terminal voltage [V]');
    title('Terminal voltage');
end

function plot_overlay_panel(traces, colors, field_name, y_label, plot_title)
    nexttile;
    hold on;
    for i = 1:numel(traces)
        tr = traces(i);
        if isfield(tr, field_name)
            plot(tr.time_hr, tr.(field_name), 'Color', colors(i, :), 'LineWidth', 0.8);
        end
    end
    if contains(y_label, 'mol/m^3')
        yline(0, ':r');
    end
    grid on;
    xlabel('Checkpoint sim time [h]');
    ylabel(y_label);
    title(plot_title);
end

function plot_terminal_panel(x, traces, field_name, plot_title)
    y = nan(numel(traces), 1);
    for i = 1:numel(traces)
        if isfield(traces(i), field_name) && ~isempty(traces(i).(field_name))
            y(i) = traces(i).(field_name)(end);
        end
    end
    plot(x, y, 'o-', 'LineWidth', 1.0);
    if contains(plot_title, 'ce') || contains(plot_title, 'log')
        yline(0, ':r');
    end
    grid on;
    xlabel('Checkpoint segment index');
    title(plot_title);
end

function add_checkpoint_colorbar(colors, n)
    cb = colorbar;
    cb.Layout.Tile = 'east';
    colormap(colors);
    if n > 1
        caxis([1, n]);
    else
        caxis([0.5, 1.5]);
    end
    cb.Label.String = 'manifest checkpoint order';
end

function plot_checkpoint_manifest_overlay(overlay, failure_trace)
    n = numel(overlay.traces);
    colors = parula(max(n, 2));
    figure('Name', 'Checkpoint manifest overlay: last-segment electrolyte trajectories', 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    hold on;
    for i = 1:n
        tr = overlay.traces(i);
        plot(tr.elapsed_hr, tr.ce0n, 'Color', colors(i, :), 'LineWidth', 0.9);
    end
    yline(0, ':r');
    grid on;
    xlabel('Elapsed in checkpoint last segment [h]');
    ylabel('ce0n [mol/m^3]');
    title('Anode boundary electrolyte at checkpoint segments');

    nexttile;
    hold on;
    for i = 1:n
        tr = overlay.traces(i);
        plot(tr.elapsed_hr, tr.log_arg_min, 'Color', colors(i, :), 'LineWidth', 0.9);
    end
    yline(0, ':r');
    grid on;
    xlabel('Elapsed in checkpoint last segment [h]');
    ylabel('min log argument [mol/m^3]');
    title('Minimum electrolyte log argument');

    nexttile;
    hold on;
    for i = 1:n
        tr = overlay.traces(i);
        plot(tr.elapsed_hr, tr.c_e_min, 'Color', colors(i, :), 'LineWidth', 0.9);
    end
    yline(0, ':r');
    grid on;
    xlabel('Elapsed in checkpoint last segment [h]');
    ylabel('min internal c_e [mol/m^3]');
    title('Internal electrolyte minimum');

    nexttile;
    hold on;
    for i = 1:n
        tr = overlay.traces(i);
        plot(tr.elapsed_hr, tr.epsilon_e_n, 'Color', colors(i, :), 'LineWidth', 0.9);
    end
    yline(0, ':r');
    grid on;
    xlabel('Elapsed in checkpoint last segment [h]');
    ylabel('\epsilon_{e,n} [-]');
    title('Anode porosity');

    cb = colorbar;
    cb.Layout.Tile = 'east';
    colormap(colors);
    if n > 1
        caxis([1, n]);
    else
        caxis([0.5, 1.5]);
    end
    cb.Label.String = 'manifest checkpoint order';

    figure('Name', 'Checkpoint manifest terminal electrolyte trend', 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    terminal_ce0n = nan(n, 1);
    terminal_log_arg = nan(n, 1);
    terminal_cemin = nan(n, 1);
    terminal_eps = nan(n, 1);
    for i = 1:n
        terminal_ce0n(i) = overlay.traces(i).ce0n(end);
        terminal_log_arg(i) = overlay.traces(i).log_arg_min(end);
        terminal_cemin(i) = overlay.traces(i).c_e_min(end);
        terminal_eps(i) = overlay.traces(i).epsilon_e_n(end);
    end
    x = overlay.segment_idx(:);

    nexttile;
    plot(x, terminal_ce0n, 'o-r', 'LineWidth', 1.0); grid on;
    yline(0, ':r');
    xlabel('Checkpoint segment index');
    ylabel('terminal ce0n');
    title('Terminal ce0n');

    nexttile;
    plot(x, terminal_log_arg, 'o-k', 'LineWidth', 1.0); grid on;
    yline(0, ':r');
    xlabel('Checkpoint segment index');
    ylabel('terminal min log arg');
    title('Terminal min log argument');

    nexttile;
    plot(x, terminal_cemin, 'o-b', 'LineWidth', 1.0); grid on;
    yline(0, ':r');
    xlabel('Checkpoint segment index');
    ylabel('terminal min internal c_e');
    title('Terminal internal electrolyte minimum');

    nexttile;
    plot(x, terminal_eps, 'o-g', 'LineWidth', 1.0); grid on;
    xlabel('Checkpoint segment index');
    ylabel('terminal \epsilon_{e,n}');
    title('Terminal anode porosity');

    if nargin >= 2 && isstruct(failure_trace)
        fprintf('\nCheckpoint overlay loaded %d trajectories from manifest:\n  %s\n', ...
            n, overlay.manifest_file);
    end
end

function manifest_path = checkpoint_manifest_path_from_checkpoint(checkpoint_file)
    [folder, base_name, ~] = fileparts(char(checkpoint_file));
    token = regexp(base_name, '^(.*)_repeat\d+_segment\d+$', 'tokens', 'once');
    if isempty(token)
        manifest_path = '';
    else
        manifest_path = fullfile(folder, [token{1} '_manifest.mat']);
    end
end

function plot_trace(trace)
    t_hr = trace.time / 3600;
    fail_t = NaN;
    if isfield(trace, 'failure') && trace.failure.occurred && ~isempty(trace.time)
        fail_t = trace.time(end) / 3600;
    end
    idx_nan_focus = first_nan_focus_idx(trace);
    if isempty(idx_nan_focus)
        nan_focus_t = fail_t;
    else
        nan_focus_t = trace.time(idx_nan_focus) / 3600;
    end

    figure('Name', 'SPMe warm-start failure trace: state consistency', 'Color', 'w');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, trace.current, 'm-', 'LineWidth', 1.0); grid on;
    ylabel('Current [A]');
    title('Current');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.voltage, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('Voltage [V]');
    title('Voltage');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.epsilon_e_n_output, 'g-', 'LineWidth', 1.2); hold on;
    plot(t_hr, trace.epsilon_e_n_state, 'k--', 'LineWidth', 1.0);
    yline(0, ':r');
    grid on;
    ylabel('\epsilon_{e,n} [-]');
    legend('optional output', 'raw state', 'Location', 'best');
    title('Anode Porosity Consistency');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.c_e_min, 'r-', 'LineWidth', 1.1); hold on;
    plot(t_hr, trace.c_e_mean, 'k-', 'LineWidth', 1.0);
    plot(t_hr, trace.c_e_max, 'b-', 'LineWidth', 1.0);
    yline(0, ':r');
    grid on;
    ylabel('c_e [mol/m^3]');
    legend('min', 'mean', 'max', 'Location', 'best');
    title('Electrolyte Concentration');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.delta_sei * 1e9, 'r-', 'LineWidth', 1.0); hold on;
    plot(t_hr, trace.delta_li * 1e9, 'b-', 'LineWidth', 1.0);
    plot(t_hr, trace.delta_film * 1e9, 'k--', 'LineWidth', 1.0);
    grid on;
    ylabel('Thickness [nm]');
    legend('\delta_{sei}', '\delta_{Li}', '\delta_{film}', 'Location', 'best');
    title('Film States');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.c_ec_sei_min, 'r-', 'LineWidth', 1.1); hold on;
    plot(t_hr, trace.c_ec_sei_mean, 'k-', 'LineWidth', 1.0);
    plot(t_hr, trace.c_ec_sei_max, 'b-', 'LineWidth', 1.0);
    yline(0, ':r');
    grid on;
    ylabel('c_{EC,SEI} [mol/m^3]');
    legend('min', 'mean', 'max', 'Location', 'best');
    title('Safari EC States');
    mark_failure(fail_t);

    figure('Name', 'SPMe warm-start failure trace: electrolyte outputs', 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, trace.soc, 'b-', 'LineWidth', 1.0); grid on;
    ylabel('SoC [-]');
    title('SoC');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.R_film, 'Color', [0.5 0 0.5], 'LineWidth', 1.0); grid on;
    ylabel('R_{film} [\Omega m^2]');
    title('Film Resistance');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.V_electrolytePolar, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('V_{electrolytePolar} [V]');
    title('Electrolyte Polarization');
    mark_failure(fail_t);

    nexttile;
    plot(t_hr, trace.cen_bar, 'r-', 'LineWidth', 1.0); hold on;
    plot(t_hr, trace.ces_bar, 'k-', 'LineWidth', 1.0);
    plot(t_hr, trace.cep_bar, 'b-', 'LineWidth', 1.0);
    yline(0, ':r');
    grid on;
    ylabel('average c_e [mol/m^3]');
    legend('anode', 'separator', 'cathode', 'Location', 'best');
    title('Regional Average Electrolyte');
    mark_failure(fail_t);

    figure('Name', 'SPMe warm-start failure trace: electrolyte log arguments near NaN onset', 'Color', 'w');
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, trace.ce0n, 'r-', 'LineWidth', 1.0); hold on;
    plot(t_hr, trace.cens, 'Color', [0.85 0.45 0], 'LineWidth', 1.0);
    plot(t_hr, trace.cesp, 'g-', 'LineWidth', 1.0);
    plot(t_hr, trace.ce0p, 'b-', 'LineWidth', 1.0);
    yline(0, ':r');
    grid on;
    ylabel('boundary c_e [mol/m^3]');
    legend('ce0n', 'cens', 'cesp', 'ce0p', 'Location', 'best');
    title('Arguments to log() in V_{electrolytePolar}');
    mark_nan_focus(nan_focus_t);
    zoom_to_focus(nan_focus_t, t_hr);

    nexttile;
    plot(t_hr, trace.electrolyte_log_arg_min, 'k-', 'LineWidth', 1.1); hold on;
    yline(0, ':r');
    grid on;
    xlabel('Replay time [h]');
    ylabel('min log argument [mol/m^3]');
    title('First crossing below zero explains NaN electrolyte polarization');
    mark_nan_focus(nan_focus_t);
    zoom_to_focus(nan_focus_t, t_hr);

    figure('Name', 'SPMe warm-start failure trace: voltage NaN onset zoom', 'Color', 'w');
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_hr, trace.voltage, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('Voltage [V]');
    title('Terminal Voltage Near First NaN');
    mark_nan_focus(nan_focus_t);
    zoom_to_focus(nan_focus_t, t_hr);

    nexttile;
    plot(t_hr, trace.V_electrolytePolar, 'k-', 'LineWidth', 1.0); grid on;
    ylabel('V_{electrolytePolar} [V]');
    title('Electrolyte Polarization Near First NaN');
    mark_nan_focus(nan_focus_t);
    zoom_to_focus(nan_focus_t, t_hr);

    nexttile;
    plot(t_hr, trace.electrolyte_log_arg_min, 'k-', 'LineWidth', 1.1); hold on;
    yline(0, ':r');
    grid on;
    xlabel('Replay time [h]');
    ylabel('min log argument [mol/m^3]');
    title('Minimum Electrolyte Log Argument Near First NaN');
    mark_nan_focus(nan_focus_t);
    zoom_to_focus(nan_focus_t, t_hr);
end

function idx = first_nan_focus_idx(trace)
    candidates = [];
    idx_voltage = first_nonfinite_after_finite(trace.voltage);
    idx_vpolar = first_nonfinite_after_finite(trace.V_electrolytePolar);
    idx_log_arg = first_bad_log_arg_after_positive(trace.electrolyte_log_arg_min);
    if ~isempty(idx_voltage)
        candidates(end+1, 1) = idx_voltage; %#ok<AGROW>
    end
    if ~isempty(idx_vpolar)
        candidates(end+1, 1) = idx_vpolar; %#ok<AGROW>
    end
    if ~isempty(idx_log_arg)
        candidates(end+1, 1) = idx_log_arg; %#ok<AGROW>
    end
    if isempty(candidates)
        idx = [];
    else
        idx = min(candidates);
    end
end

function idx = first_nonfinite_after_finite(v)
    idx = [];
    v = v(:);
    first_finite = find(isfinite(v), 1, 'first');
    if isempty(first_finite)
        return;
    end
    rel = find(~isfinite(v(first_finite:end)), 1, 'first');
    if ~isempty(rel)
        idx = first_finite + rel - 1;
    end
end

function idx = first_bad_log_arg_after_positive(v)
    idx = [];
    v = v(:);
    first_positive = find(isfinite(v) & v > 0, 1, 'first');
    if isempty(first_positive)
        return;
    end
    rel = find(v(first_positive:end) <= 0 | ~isfinite(v(first_positive:end)), 1, 'first');
    if ~isempty(rel)
        idx = first_positive + rel - 1;
    end
end

function zoom_to_focus(focus_t_hr, t_hr)
    if ~isfinite(focus_t_hr) || isempty(t_hr)
        return;
    end
    half_width_hr = 0.15;
    x_min = max(min(t_hr), focus_t_hr - half_width_hr);
    x_max = min(max(t_hr), focus_t_hr + half_width_hr);
    if x_max > x_min
        xlim([x_min, x_max]);
    end
end

function mark_nan_focus(focus_t_hr)
    if isfinite(focus_t_hr)
        xline(focus_t_hr, '--r', 'first NaN/invalid log', 'HandleVisibility', 'off');
    end
end

function print_nan_origin_summary(trace)
    fprintf('\nWarm-start NaN origin summary:\n');
    idx_vnan = find(~isfinite(trace.voltage), 1, 'first');
    idx_vpol_nan = find(~isfinite(trace.V_electrolytePolar), 1, 'first');
    idx_log_bad = first_bad_log_arg_after_positive(trace.electrolyte_log_arg_min);
    print_event('First nonfinite terminal voltage', trace, idx_vnan);
    print_event('First nonfinite V_electrolytePolar', trace, idx_vpol_nan);
    print_event('First invalid electrolyte log argument', trace, idx_log_bad);
    if ~isempty(idx_log_bad)
        fprintf('  ce0n=%.12g, cens=%.12g, cesp=%.12g, ce0p=%.12g\n', ...
            trace.ce0n(idx_log_bad), trace.cens(idx_log_bad), ...
            trace.cesp(idx_log_bad), trace.ce0p(idx_log_bad));
        [~, which_min] = min([trace.ce0n(idx_log_bad), trace.cens(idx_log_bad), ...
            trace.cesp(idx_log_bad), trace.ce0p(idx_log_bad)]);
        names = {'ce0n', 'cens', 'cesp', 'ce0p'};
        fprintf('  first bad log argument appears at %s.\n', names{which_min});
    end
    fprintf('\n');
end

function print_event(label, trace, idx)
    if isempty(idx)
        fprintf('  %s: not found\n', label);
        return;
    end
    fprintf('  %s: sample %d, replay t=%.6g h, original segment=%d, local segment=%d, V=%.12g, V_electrolytePolar=%.12g\n', ...
        label, idx, trace.time(idx) / 3600, trace.original_segment_idx(idx), ...
        trace.local_segment_idx(idx), trace.voltage(idx), trace.V_electrolytePolar(idx));
end

function mark_failure(fail_t)
    if isfinite(fail_t)
        xline(fail_t, '--r', 'failure', 'HandleVisibility', 'off');
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

function x = get_field_or(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        x = s.(name);
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

function ensure_test_paths()
    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));
end

function req_outputs = default_req_outputs()
    req_outputs = {'delta_sei', 'delta_sei_dot', ...
        'delta_li', 'delta_li_dot', 'delta_film', 'delta_film_dot', ...
        'epsilon_e_n', 'epsilon_e_n_dot', ...
        'j_sei', 'j_plating', 'j_side', 'R_film', ...
        'V_electrolytePolar', 'ce0n', 'cens', 'cesp', 'ce0p', ...
        'cen_bar', 'ces_bar', 'cep_bar'};
end

function tf = is_manifest_file(path_value)
    tf = false;
    if ~(ischar(path_value) || isstring(path_value))
        return;
    end
    path_value = char(path_value);
    if ~isfile(path_value)
        return;
    end
    [~, name, ext] = fileparts(path_value);
    tf = strcmpi(ext, '.mat') && endsWith(name, '_manifest');
end

function v = last_or_nan(x)
    if isempty(x)
        v = NaN;
    else
        v = x(end);
    end
end

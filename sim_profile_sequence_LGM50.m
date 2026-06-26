function sim = sim_profile_sequence_LGM50(model_name, init_soc, profile_sequence, varargin)
%SIM_PROFILE_SEQUENCE_LGM50 Run concatenated profile segments continuously.
%   sim = sim_profile_sequence_LGM50(model_name, init_soc, profile_sequence, ...)
%
%   profile_sequence examples:
%     'FUDS'
%     {'CC', -1}
%     {'CC', -1, 'CC', 1}
%     {'FUDS', 'CC', -1}
%     {'REST', 0.1, 'CC', -1}   % REST duration is in hours
%
%   Name-Value options:
%     'Repeat' : positive integer, repeat full sequence N times (default 1)
%     'VMin'   : custom lower cutoff voltage [V] (default [], disabled)
%     'VMax'   : custom upper cutoff voltage [V] (default [], disabled)
%     'ShowProgress' : true/false for waitbar display (default true)
%     'MemoryLight' : true/false to avoid storing large per-step state
%                     histories (x/y/c_e/etc.) in the aggregated output
%                     sequence; keeps only core signals + requested
%                     optional outputs (default false)
%     'InitialState' : state struct for starting the first segment from a
%                     previously saved terminal state
%     'CheckpointFile' : optional .mat path saved after each completed
%                     repeat with aggregate output and warm-start state
%     'EnableDiagnostics' : true/false to save segment handoff/failure
%                     diagnostics (default false)
%     'DiagnosticsFile' : optional .mat path for sequence diagnostics
%     'DiagnosticContext' : optional caller metadata copied into diagnostics
%     'SPMeParamOverrides' : optional struct passed to sim_profile_LGM50

    if nargin < 1 || isempty(model_name)
        model_name = 'SPMe';
    end
    if nargin < 2 || isempty(init_soc)
        init_soc = 0.9;
    end
    if nargin < 3 || isempty(profile_sequence)
        profile_sequence = 'DST';
    end

    parser = inputParser;
    addParameter(parser, 'Nxn',10);
    addParameter(parser, 'Nxs',5);
    addParameter(parser, 'Nxp',10);
    addParameter(parser, 'Nr', 15);
    addParameter(parser, 'DeltaT', 1);
    addParameter(parser, 'Tamb', 25);
    addParameter(parser, 'SaveOutput', false);
    addParameter(parser, 'OutputFile', '');
    addParameter(parser, 'ECMOrder', 2);
    addParameter(parser, 'Crate', 1);
    addParameter(parser, 'DegradationModel', 'none');
    addParameter(parser, 'ElectrolyteModel', 'fdm');
    addParameter(parser, 'EnableLithiumPlating', []);
    addParameter(parser, 'EnableAnodePorosityLoss', []);
    addParameter(parser, 'DFNLplCouplingMode', 'spme_compatible');
    addParameter(parser, 'VEnd', []);
    addParameter(parser, 'Repeat', 1);
    addParameter(parser, 'VMin', []);
    addParameter(parser, 'VMax', []);
    addParameter(parser, 'ShowProgress', true);
    addParameter(parser, 'SPMeOutputs', {});
    addParameter(parser, 'MemoryLight', false);
    addParameter(parser, 'InitialState', struct());
    addParameter(parser, 'CheckpointFile', '');
    addParameter(parser, 'EnableDiagnostics', false);
    addParameter(parser, 'DiagnosticsFile', '');
    addParameter(parser, 'DiagnosticContext', struct());
    addParameter(parser, 'SPMeParamOverrides', struct());
    parse(parser, varargin{:});
    opt = parser.Results;
    if ~isscalar(opt.Repeat) || ~isfinite(opt.Repeat) || opt.Repeat < 1 || abs(opt.Repeat - round(opt.Repeat)) > 0
        error('Repeat must be a positive integer.');
    end
    if ~isempty(opt.VMin) && (~isscalar(opt.VMin) || ~isfinite(opt.VMin))
        error('VMin must be a finite scalar when provided.');
    end
    if ~isempty(opt.VMax) && (~isscalar(opt.VMax) || ~isfinite(opt.VMax))
        error('VMax must be a finite scalar when provided.');
    end
    if ~isempty(opt.VMin) && ~isempty(opt.VMax) && opt.VMin >= opt.VMax
        error('VMin must be smaller than VMax.');
    end
    opt.EnableDiagnostics = logical(opt.EnableDiagnostics);

    segments = parse_profile_sequence(profile_sequence);
    segments = repeat_segments(segments, opt.Repeat);
    model_name = upper(char(model_name));
    if ~ismember(model_name, {'DFN', 'SPME', 'ECM'})
        error('model_name must be ''DFN'', ''SPMe'', or ''ECM''.');
    end
    if opt.EnableDiagnostics && isempty(opt.DiagnosticsFile)
        opt.DiagnosticsFile = default_sequence_diagnostics_file(opt.CheckpointFile, model_name);
    end

    n_segments_total = numel(segments);
    capacity_cycle_by_segment = assign_capacity_cycle_index(segments);
    h_wait = [];
    if opt.ShowProgress
        try
            h_wait = waitbar(0, sprintf('Running %s sequence (0/%d segments)...', model_name, n_segments_total), ...
                'Name', 'sim_profile_sequence_LGM50');
        catch
            h_wait = [];
        end
    end
    cleanup_wait = onCleanup(@() close_waitbar_safe(h_wait)); %#ok<NASGU>

    initial_state = struct();
    if isstruct(opt.InitialState) && ~isempty(fieldnames(opt.InitialState))
        initial_state = opt.InitialState;
    end
    save_sequence_diagnostic_event(opt, 'sequence_start', model_name, init_soc, ...
        1, n_segments_total, segments, initial_state, []);
    seg0 = run_segment(model_name, init_soc, segments(1), opt, initial_state, 1, n_segments_total);
    seg0.sim = apply_voltage_window(seg0.sim, opt.VMin, opt.VMax);
    save_sequence_diagnostic_event(opt, 'segment_complete', model_name, init_soc, ...
        1, n_segments_total, segments, initial_state, seg0.sim);
    update_waitbar_safe(h_wait, 1, n_segments_total, segments(1));
    sim = seg0.sim;
    if opt.MemoryLight
        sim = prune_heavy_history_fields(sim);
    end
    sim.profile = 'SEQUENCE';
    sim.profile_definition = segments;
    sim.profile_sequence_input = profile_sequence;
    sim.repeat_count = opt.Repeat;
    sim.segment_index = ones(numel(sim.time), 1);
    sim.segment_cycle = ones(numel(sim.time), 1);
    sim.capacity_cycle = capacity_cycle_by_segment(1) * ones(numel(sim.time), 1);
    [charged_ah0, discharged_ah0] = compute_charge_discharge_ah(seg0.sim.time, seg0.sim.current);
    sim.segment_log = struct('segment', 1, ...
        'cycle', segments(1).cycle, ...
        'capacity_cycle', capacity_cycle_by_segment(1), ...
        'segment_in_cycle', segments(1).segment_in_cycle, ...
        'name', segments(1).name, ...
        'crate', segments(1).crate, 'terminated_early', seg0.sim.terminated_early, ...
        'termination_reason', seg0.sim.termination_reason, ...
        'charged_Ah', charged_ah0, ...
        'discharged_Ah', discharged_ah0, ...
        'net_Ah', discharged_ah0 - charged_ah0);
    if is_end_of_repeat(1, segments)
        save_sequence_checkpoint(opt.CheckpointFile, sim, seg0.sim, model_name, ...
            1, n_segments_total, segments, opt);
    end

    for s = 2:numel(segments)
        init_state = terminal_state_from_sim(seg0.sim, model_name);
        init_soc_seg = seg0.sim.soc(end);
        save_sequence_diagnostic_event(opt, 'segment_start', model_name, init_soc_seg, ...
            s, n_segments_total, segments, init_state, []);
        seg0 = run_segment(model_name, init_soc_seg, segments(s), opt, init_state, s, n_segments_total);
        seg0.sim = apply_voltage_window(seg0.sim, opt.VMin, opt.VMax);
        save_sequence_diagnostic_event(opt, 'segment_complete', model_name, init_soc_seg, ...
            s, n_segments_total, segments, init_state, seg0.sim);
        update_waitbar_safe(h_wait, s, n_segments_total, segments(s));
        seg = seg0.sim;
        idx = 2:numel(seg.time);
        if isempty(idx)
            continue;
        end

        t_offset = sim.time(end);
        sim.time = [sim.time; t_offset + seg.time(idx)];
        sim.current = [sim.current; seg.current(idx)];
        sim.voltage = [sim.voltage; seg.voltage(idx)];
        sim.soc = [sim.soc; seg.soc(idx)];
        sim.profile_time = [sim.profile_time; seg.profile_time(idx)];
        sim.segment_index = [sim.segment_index; s*ones(numel(idx),1)];
        sim.segment_cycle = [sim.segment_cycle; segments(s).cycle*ones(numel(idx),1)];
        sim.capacity_cycle = [sim.capacity_cycle; capacity_cycle_by_segment(s)*ones(numel(idx),1)];

        sim = append_optional_fields(sim, seg, idx, opt.MemoryLight);
        [charged_ah_seg, discharged_ah_seg] = compute_charge_discharge_ah(seg.time, seg.current);
        sim.segment_log(end+1) = struct('segment', s, ...
            'cycle', segments(s).cycle, ...
            'capacity_cycle', capacity_cycle_by_segment(s), ...
            'segment_in_cycle', segments(s).segment_in_cycle, ...
            'name', segments(s).name, ...
            'crate', segments(s).crate, 'terminated_early', seg.terminated_early, ...
            'termination_reason', seg.termination_reason, ...
            'charged_Ah', charged_ah_seg, ...
            'discharged_Ah', discharged_ah_seg, ...
            'net_Ah', discharged_ah_seg - charged_ah_seg); %#ok<AGROW>
        sim.t_last = numel(sim.time);
        if is_end_of_repeat(s, segments)
            save_sequence_checkpoint(opt.CheckpointFile, sim, seg0.sim, model_name, ...
                s, n_segments_total, segments, opt);
        end
    end

    sim.terminated_early = any([sim.segment_log.terminated_early]);
    if sim.terminated_early
        sim.termination_reason = 'At least one segment terminated early by voltage condition.';
    else
        sim.termination_reason = 'Sequence completed.';
    end
    if ~isempty(opt.VMin) || ~isempty(opt.VMax)
        sim.voltage_window = struct('VMin', opt.VMin, 'VMax', opt.VMax);
    end
    sim.capacity_cycle_log = build_capacity_cycle_log(sim.segment_log);
    sim.saved_file = '';

    if opt.SaveOutput
        base_dir = fileparts(mfilename('fullpath'));
        if isempty(opt.OutputFile)
            out_dir = fullfile(base_dir, 'results', model_name);
            if ~isfolder(out_dir)
                mkdir(out_dir);
            end
            out_name = sprintf('%s_SEQUENCE_%0.3fSOC.mat', lower(model_name), init_soc);
            out_path = fullfile(out_dir, out_name);
        else
            out_path = opt.OutputFile;
        end
        save(out_path, 'sim');
        sim.saved_file = out_path;
    end
    save_sequence_diagnostic_event(opt, 'sequence_complete', model_name, init_soc, ...
        n_segments_total, n_segments_total, segments, terminal_state_from_sim(seg0.sim, model_name), seg0.sim);
end

function segments_rep = repeat_segments(segments_base, repeat_count)
    segments_rep = struct('name', {}, 'crate', {}, 'rest_hours', {}, 'cycle', {}, 'segment_in_cycle', {});
    n = numel(segments_base);
    for c = 1:repeat_count
        for s = 1:n
            segments_rep(end+1) = struct( ... %#ok<AGROW>
                'name', segments_base(s).name, ...
                'crate', segments_base(s).crate, ...
                'rest_hours', segments_base(s).rest_hours, ...
                'cycle', c, ...
                'segment_in_cycle', s);
        end
    end
end

function cap_cycle = assign_capacity_cycle_index(segments)
    n = numel(segments);
    cap_cycle = ones(n, 1);
    if n <= 1
        return;
    end

    seg_dir = zeros(n, 1);
    for i = 1:n
        if strcmp(segments(i).name, 'CC')
            seg_dir(i) = sign(segments(i).crate);
        elseif strcmp(segments(i).name, 'REST')
            seg_dir(i) = 0;
        else
            % Dynamic drive cycles are generally discharge-oriented.
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

function [charged_ah, discharged_ah] = compute_charge_discharge_ah(time_s, current_a)
    if isempty(time_s) || numel(time_s) < 2 || isempty(current_a)
        charged_ah = 0;
        discharged_ah = 0;
        return;
    end
    t = time_s(:);
    i = current_a(:);
    i_charge = max(-i, 0);
    i_discharge = max(i, 0);
    charged_ah = trapz(t, i_charge) / 3600;
    discharged_ah = trapz(t, i_discharge) / 3600;
end

function cycle_log = build_capacity_cycle_log(segment_log)
    cycle_ids = [segment_log.capacity_cycle];
    unique_cycles = unique(cycle_ids, 'stable');
    cycle_log = struct('cycle', {}, 'start_segment', {}, 'end_segment', {}, ...
        'charged_Ah', {}, 'discharged_Ah', {}, 'net_Ah', {});

    for k = 1:numel(unique_cycles)
        cyc = unique_cycles(k);
        idx = find(cycle_ids == cyc);
        charged_ah = sum([segment_log(idx).charged_Ah]);
        discharged_ah = sum([segment_log(idx).discharged_Ah]);
        cycle_log(end+1) = struct( ... %#ok<AGROW>
            'cycle', cyc, ...
            'start_segment', segment_log(idx(1)).segment, ...
            'end_segment', segment_log(idx(end)).segment, ...
            'charged_Ah', charged_ah, ...
            'discharged_Ah', discharged_ah, ...
            'net_Ah', discharged_ah - charged_ah);
    end
end

function out = run_segment(model_name, init_soc, seg, opt, init_state, segment_idx, n_segments_total)
    if strcmp(seg.name, 'CC')
        prof = {'CC', seg.crate};
    elseif strcmp(seg.name, 'REST')
        prof = {'REST', seg.rest_hours};
    else
        prof = seg.name;
    end
    tmp_sim = [];
    diagnostic_context = segment_diagnostic_context(model_name, init_soc, seg, opt, init_state, segment_idx, n_segments_total);
    step_diagnostics_file = segment_step_diagnostics_file(opt, seg, segment_idx);
    try
        evalc(['tmp_sim = sim_profile_LGM50(model_name, init_soc, prof, ' ...
            '''Nxn'', opt.Nxn, ''Nxs'', opt.Nxs, ''Nxp'', opt.Nxp, ''Nr'', opt.Nr, ' ...
            '''DeltaT'', opt.DeltaT, ''Tamb'', opt.Tamb, ''ECMOrder'', opt.ECMOrder, ' ...
            '''Crate'', opt.Crate, ''DegradationModel'', opt.DegradationModel, ' ...
            '''ElectrolyteModel'', opt.ElectrolyteModel, ' ...
            '''EnableLithiumPlating'', opt.EnableLithiumPlating, ' ...
            '''EnableAnodePorosityLoss'', opt.EnableAnodePorosityLoss, ' ...
            '''DFNLplCouplingMode'', opt.DFNLplCouplingMode, ' ...
            '''VEnd'', opt.VEnd, ''InitialState'', init_state, ' ...            '''SuppressSummary'', true, ''SPMeOutputs'', opt.SPMeOutputs, ''SaveOutput'', false, ' ...
            '''DiagnosticsFile'', step_diagnostics_file, ''DiagnosticContext'', diagnostic_context, ' ...
            '''SPMeParamOverrides'', opt.SPMeParamOverrides);']);
    catch ME
        save_sequence_diagnostic_failure(opt, model_name, init_soc, seg, init_state, ...
            segment_idx, n_segments_total, ME, step_diagnostics_file);
        rethrow(ME);
    end
    out.sim = tmp_sim;
end

function segments = parse_profile_sequence(profile_sequence)
    valid_dynamic = {'BJDST', 'DST', 'FUDS', 'US06'};
    if ischar(profile_sequence) || isstring(profile_sequence)
        name = upper(char(profile_sequence));
        if strcmp(name, 'CC')
            error('For sequence mode, CC must be specified as {''CC'', crate}.');
        end
        if strcmp(name, 'REST')
            error('For sequence mode, REST must be specified as {''REST'', duration_hours}.');
        end
        if ~ismember(name, valid_dynamic)
            error('Unknown profile token: %s', name);
        end
        segments = struct('name', name, 'crate', NaN, 'rest_hours', NaN);
        return;
    end

    if ~iscell(profile_sequence)
        error('profile_sequence must be char/string or a cell token sequence.');
    end

    segments = struct('name', {}, 'crate', {}, 'rest_hours', {});
    i = 1;
    while i <= numel(profile_sequence)
        tok = profile_sequence{i};
        if ~(ischar(tok) || isstring(tok))
            error('Expected profile name token at position %d.', i);
        end
        name = upper(char(tok));
        if strcmp(name, 'CC')
            if i+1 > numel(profile_sequence) || ~isnumeric(profile_sequence{i+1}) || ...
                    ~isscalar(profile_sequence{i+1}) || profile_sequence{i+1} == 0
                error('CC token at position %d must be followed by non-zero scalar Crate.', i);
            end
            segments(end+1) = struct('name', 'CC', 'crate', profile_sequence{i+1}, 'rest_hours', NaN); %#ok<AGROW>
            i = i + 2;
        elseif strcmp(name, 'REST')
            if i+1 > numel(profile_sequence) || ~isnumeric(profile_sequence{i+1}) || ...
                    ~isscalar(profile_sequence{i+1}) || ~isfinite(profile_sequence{i+1}) || profile_sequence{i+1} <= 0
                error('REST token at position %d must be followed by positive scalar duration in hours.', i);
            end
            segments(end+1) = struct('name', 'REST', 'crate', 0, 'rest_hours', profile_sequence{i+1}); %#ok<AGROW>
            i = i + 2;
        elseif ismember(name, valid_dynamic)
            segments(end+1) = struct('name', name, 'crate', NaN, 'rest_hours', NaN); %#ok<AGROW>
            i = i + 1;
        else
            error('Unknown profile token at position %d: %s', i, name);
        end
    end
    if isempty(segments)
        error('profile_sequence cannot be empty.');
    end
end

function init_state = terminal_state_from_sim(sim, model_name)
    init_state = struct();
    switch upper(char(model_name))
        case 'SPME'
            init_state.x = sim.x(:, end);
            init_state.p = sim.p;
        case 'DFN'
            init_state.x = sim.x(:, end);
            init_state.z = sim.z(:, end);
            init_state.p = sim.p;
        case 'ECM'
            init_state.x = sim.x(:, end);
            init_state.p = sim.p;
    end
end

function tf = is_end_of_repeat(segment_idx, segments)
    tf = segment_idx >= numel(segments) || segments(segment_idx + 1).cycle ~= segments(segment_idx).cycle;
end

function save_sequence_checkpoint(checkpoint_file, sim, last_segment, model_name, segment_idx, n_segments_total, segments, opt)
    if isempty(checkpoint_file)
        return;
    end
    checkpoint_file = char(checkpoint_file);
    checkpoint_path = checkpoint_path_for_segment(checkpoint_file, segment_idx, segments);
    [checkpoint_dir, ~, ~] = fileparts(checkpoint_path);
    if ~isempty(checkpoint_dir) && ~isfolder(checkpoint_dir)
        mkdir(checkpoint_dir);
    end

    checkpoint = struct();
    checkpoint.saved_at = datestr(now, 31);
    checkpoint.model_name = model_name;
    checkpoint.segment_idx = segment_idx;
    checkpoint.n_segments_total = n_segments_total;
    checkpoint.next_segment_idx = min(segment_idx + 1, n_segments_total + 1);
    checkpoint.completed_fraction = segment_idx / max(n_segments_total, 1);
    checkpoint.checkpoint_file = checkpoint_path;
    checkpoint.checkpoint_manifest_file = checkpoint_manifest_path(checkpoint_file);
    checkpoint.segments = segments;
    checkpoint.options = opt;
    checkpoint.sim = sim;
    checkpoint.last_segment = last_segment;
    checkpoint.terminal_state = terminal_state_from_sim(last_segment, model_name);
    checkpoint.terminal_soc = last_segment.soc(end);
    checkpoint.terminal_voltage = last_segment.voltage(end);
    checkpoint.terminal_time = sim.time(end);
    if isfield(sim, 'capacity_cycle_log')
        checkpoint.capacity_cycle_log = sim.capacity_cycle_log;
    else
        checkpoint.capacity_cycle_log = build_capacity_cycle_log(sim.segment_log);
    end
    save(checkpoint_path, 'checkpoint', '-v7.3');
    update_checkpoint_manifest(checkpoint_file, checkpoint_path, checkpoint);
end

function checkpoint_path = checkpoint_path_for_segment(checkpoint_file, segment_idx, segments)
    [folder, base_name, ext] = fileparts(checkpoint_file);
    if isempty(ext)
        ext = '.mat';
    end
    repeat_idx = segments(segment_idx).cycle;
    checkpoint_name = sprintf('%s_repeat%04d_segment%05d%s', ...
        base_name, repeat_idx, segment_idx, ext);
    checkpoint_path = fullfile(folder, checkpoint_name);
end

function manifest_path = checkpoint_manifest_path(checkpoint_file)
    [folder, base_name, ~] = fileparts(checkpoint_file);
    manifest_path = fullfile(folder, [base_name '_manifest.mat']);
end

function update_checkpoint_manifest(checkpoint_file, checkpoint_path, checkpoint)
    manifest_path = checkpoint_manifest_path(checkpoint_file);
    checkpoint_manifest = struct();
    if isfile(manifest_path)
        loaded = load(manifest_path, 'checkpoint_manifest');
        if isfield(loaded, 'checkpoint_manifest')
            checkpoint_manifest = loaded.checkpoint_manifest;
        end
    end
    if ~isfield(checkpoint_manifest, 'base_file')
        checkpoint_manifest.base_file = checkpoint_file;
        checkpoint_manifest.files = {};
        checkpoint_manifest.segment_idx = [];
        checkpoint_manifest.repeat_idx = [];
        checkpoint_manifest.saved_at = {};
        checkpoint_manifest.terminal_soc = [];
        checkpoint_manifest.terminal_voltage = [];
        checkpoint_manifest.terminal_time = [];
    end
    checkpoint_manifest.files{end+1, 1} = checkpoint_path;
    checkpoint_manifest.segment_idx(end+1, 1) = checkpoint.segment_idx;
    checkpoint_manifest.repeat_idx(end+1, 1) = checkpoint.segments(checkpoint.segment_idx).cycle;
    checkpoint_manifest.saved_at{end+1, 1} = checkpoint.saved_at;
    checkpoint_manifest.terminal_soc(end+1, 1) = checkpoint.terminal_soc;
    checkpoint_manifest.terminal_voltage(end+1, 1) = checkpoint.terminal_voltage;
    checkpoint_manifest.terminal_time(end+1, 1) = checkpoint.terminal_time;
    save(manifest_path, 'checkpoint_manifest', '-v7.3');
end

function diagnostics_file = default_sequence_diagnostics_file(checkpoint_file, model_name)
    run_ts = datestr(now, 'yyyymmdd_HHMMSS');
    if ~isempty(checkpoint_file)
        [folder, base_name, ~] = fileparts(char(checkpoint_file));
        diagnostics_file = fullfile(folder, sprintf('%s_diagnostics_%s.mat', base_name, run_ts));
        return;
    end
    base_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(base_dir, 'results', model_name);
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    diagnostics_file = fullfile(out_dir, sprintf('%s_sequence_diagnostics_%s.mat', lower(model_name), run_ts));
end

function save_sequence_diagnostic_event(opt, event_name, model_name, init_soc, ...
        segment_idx, n_segments_total, segments, init_state, seg_sim)
    if ~isfield(opt, 'EnableDiagnostics') || ~opt.EnableDiagnostics
        return;
    end
    diagnostics_file = char(opt.DiagnosticsFile);
    if isempty(diagnostics_file)
        return;
    end

    event = struct();
    event.saved_at = datestr(now, 31);
    event.event = event_name;
    event.model_name = model_name;
    event.local_segment_idx = segment_idx;
    event.n_segments_total = n_segments_total;
    event.init_soc = init_soc;
    event.caller_context = opt.DiagnosticContext;
    if segment_idx >= 1 && segment_idx <= numel(segments)
        event.segment = segments(segment_idx);
    else
        event.segment = struct();
    end
    event.init_state = summarize_initial_state(init_state);
    if nargin >= 9 && isstruct(seg_sim) && ~isempty(fieldnames(seg_sim))
        event.segment_result = summarize_segment_result(seg_sim);
    else
        event.segment_result = struct();
    end

    append_sequence_diagnostic_event(diagnostics_file, event);
end

function save_sequence_diagnostic_failure(opt, model_name, init_soc, seg, init_state, ...
        segment_idx, n_segments_total, ME, step_diagnostics_file)
    if ~isfield(opt, 'EnableDiagnostics') || ~opt.EnableDiagnostics
        return;
    end
    diagnostics_file = char(opt.DiagnosticsFile);
    if isempty(diagnostics_file)
        return;
    end
    event = struct();
    event.saved_at = datestr(now, 31);
    event.event = 'segment_failure';
    event.model_name = model_name;
    event.local_segment_idx = segment_idx;
    event.n_segments_total = n_segments_total;
    event.init_soc = init_soc;
    event.caller_context = opt.DiagnosticContext;
    event.segment = seg;
    event.init_state = summarize_initial_state(init_state);
    event.step_diagnostics_file = step_diagnostics_file;
    event.error = exception_to_struct(ME);
    append_sequence_diagnostic_event(diagnostics_file, event);
    fprintf(2, 'Sequence diagnostics saved for failed segment:\n  %s\n', diagnostics_file);
    if ~isempty(step_diagnostics_file)
        fprintf(2, 'Step diagnostics target:\n  %s\n', step_diagnostics_file);
    end
end

function append_sequence_diagnostic_event(diagnostics_file, event)
    [diagnostics_dir, ~, ~] = fileparts(diagnostics_file);
    if ~isempty(diagnostics_dir) && ~isfolder(diagnostics_dir)
        mkdir(diagnostics_dir);
    end
    sequence_diagnostic = struct();
    if isfile(diagnostics_file)
        loaded = load(diagnostics_file, 'sequence_diagnostic');
        if isfield(loaded, 'sequence_diagnostic')
            sequence_diagnostic = loaded.sequence_diagnostic;
        end
    end
    if ~isfield(sequence_diagnostic, 'events')
        sequence_diagnostic.created_at = datestr(now, 31);
        sequence_diagnostic.events = {};
    elseif ~iscell(sequence_diagnostic.events)
        sequence_diagnostic.events = num2cell(sequence_diagnostic.events(:));
    end
    sequence_diagnostic.updated_at = datestr(now, 31);
    sequence_diagnostic.events{end+1, 1} = event;
    save(diagnostics_file, 'sequence_diagnostic', '-v7.3');
end

function context = segment_diagnostic_context(model_name, init_soc, seg, opt, init_state, ...
        segment_idx, n_segments_total)
    context = struct();
    context.model_name = model_name;
    context.init_soc = init_soc;
    context.local_segment_idx = segment_idx;
    context.n_segments_total = n_segments_total;
    if isfield(opt.DiagnosticContext, 'resume_next_segment_idx')
        context.original_segment_idx_estimate = opt.DiagnosticContext.resume_next_segment_idx + segment_idx - 1;
    end
    context.segment = seg;
    context.caller_context = opt.DiagnosticContext;
    context.sequence_diagnostics_file = opt.DiagnosticsFile;
    context.options = rmfield_if_present(opt, {'InitialState'});
    context.init_state = summarize_initial_state(init_state);
end

function diagnostics_file = segment_step_diagnostics_file(opt, seg, segment_idx)
    diagnostics_file = '';
    if ~isfield(opt, 'EnableDiagnostics') || ~opt.EnableDiagnostics || isempty(opt.DiagnosticsFile)
        return;
    end
    [folder, base_name, ~] = fileparts(char(opt.DiagnosticsFile));
    if strcmp(seg.name, 'CC')
        seg_desc = sprintf('cc_%+.3gC', seg.crate);
    elseif strcmp(seg.name, 'REST')
        seg_desc = sprintf('rest_%.3gh', seg.rest_hours);
    else
        seg_desc = lower(char(seg.name));
    end
    seg_desc = regexprep(seg_desc, '[^A-Za-z0-9_.+-]', '_');
    diagnostics_file = fullfile(folder, sprintf('%s_step_local%05d_repeat%04d_segment%03d_%s.mat', ...
        base_name, segment_idx, seg.cycle, seg.segment_in_cycle, seg_desc));
end

function s = summarize_initial_state(init_state)
    s = struct();
    if ~isstruct(init_state) || isempty(fieldnames(init_state))
        return;
    end
    if isfield(init_state, 'x') && ~isempty(init_state.x)
        s.x = summarize_vector(init_state.x);
        s.x_nonfinite_idx = find(~isfinite(init_state.x(:)));
    end
    if isfield(init_state, 'z') && ~isempty(init_state.z)
        s.z = summarize_vector(init_state.z);
        s.z_nonfinite_idx = find(~isfinite(init_state.z(:)));
    end
    if isfield(init_state, 'p') && ~isempty(init_state.p)
        s.p_summary = summarize_param_state(init_state.p, init_state);
    end
end

function s = summarize_segment_result(seg_sim)
    s = struct();
    s.t_last = get_field_or(seg_sim, 't_last', NaN);
    s.start_time_s = first_or_nan(seg_sim.time);
    s.end_time_s = last_or_nan(seg_sim.time);
    s.start_current_A = first_or_nan(seg_sim.current);
    s.end_current_A = last_or_nan(seg_sim.current);
    s.start_voltage = first_or_nan(seg_sim.voltage);
    s.end_voltage = last_or_nan(seg_sim.voltage);
    s.start_soc = first_or_nan(seg_sim.soc);
    s.end_soc = last_or_nan(seg_sim.soc);
    s.terminated_early = get_field_or(seg_sim, 'terminated_early', false);
    s.termination_reason = get_field_or(seg_sim, 'termination_reason', '');
    if isfield(seg_sim, 'x') && ~isempty(seg_sim.x)
        s.x_end = summarize_vector(seg_sim.x(:, end));
        s.x_end_nonfinite_idx = find(~isfinite(seg_sim.x(:, end)));
    end
    if isfield(seg_sim, 'p')
        s.p_summary = summarize_param_state(seg_sim.p, seg_sim);
    end
end

function s = summarize_param_state(p, source)
    s = struct();
    scalar_names = {'delta_t', 'OneC', 'volt_min', 'volt_max', 'c_s_n_max', ...
        'c_s_p_max', 'c_e', 'epsilon_e_n', 'epsilon_sei', 'Nsei'};
    for i = 1:numel(scalar_names)
        name = scalar_names{i};
        if isfield(p, name) && isscalar(p.(name))
            s.(name) = p.(name);
        end
    end
    if isfield(p, 'degradation_model')
        s.degradation_model = p.degradation_model;
    end
    if isfield(p, 'electrolyte_model')
        s.electrolyte_model = p.electrolyte_model;
    end
    if isfield(p, 'enable_lithium_plating')
        s.enable_lithium_plating = p.enable_lithium_plating;
    end
    if isfield(p, 'enable_anode_porosity_loss')
        s.enable_anode_porosity_loss = p.enable_anode_porosity_loss;
    end
    if isstruct(source)
        if isfield(source, 'delta_sei') && ~isempty(source.delta_sei)
            s.delta_sei = last_or_nan(source.delta_sei);
        end
        if isfield(source, 'delta_li') && ~isempty(source.delta_li)
            s.delta_li = last_or_nan(source.delta_li);
        end
        if isfield(source, 'epsilon_e_n') && ~isempty(source.epsilon_e_n)
            s.epsilon_e_n_state = last_or_nan(source.epsilon_e_n);
        end
    end
end

function summary = summarize_vector(v)
    v = v(:);
    finite_v = v(isfinite(v));
    summary = struct('numel', numel(v), 'nonfinite_count', numel(v) - numel(finite_v));
    if isempty(finite_v)
        summary.min = NaN;
        summary.max = NaN;
        summary.mean = NaN;
        summary.first = NaN;
        summary.last = NaN;
    else
        summary.min = min(finite_v);
        summary.max = max(finite_v);
        summary.mean = mean(finite_v);
        summary.first = v(1);
        summary.last = v(end);
    end
end

function s = exception_to_struct(ME)
    s = struct('identifier', ME.identifier, 'message', ME.message, 'stack', ME.stack);
end

function s = rmfield_if_present(s, names)
    for i = 1:numel(names)
        if isfield(s, names{i})
            s = rmfield(s, names{i});
        end
    end
end

function v = get_field_or(s, name, default_value)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = default_value;
    end
end

function v = first_or_nan(x)
    if isempty(x)
        v = NaN;
    else
        v = x(1);
    end
end

function v = last_or_nan(x)
    if isempty(x)
        v = NaN;
    else
        v = x(end);
    end
end

function sim = append_optional_fields(sim, seg, idx, memory_light)
    if nargin < 4
        memory_light = false;
    end
    fields_1d = {'delta_sei', 'delta_li', 'epsilon_e_n'};
    for k = 1:numel(fields_1d)
        fn = fields_1d{k};
        if isfield(sim, fn) && isfield(seg, fn)
            sim.(fn) = [sim.(fn); seg.(fn)(idx)];
        end
    end
    fields_2d = {'x', 'y', 'z', 'c_s_n', 'c_s_p', 'c_e', 'temperature'};
    if ~memory_light
        for k = 1:numel(fields_2d)
            fn = fields_2d{k};
            if isfield(sim, fn) && isfield(seg, fn)
                sim.(fn) = [sim.(fn), seg.(fn)(:, idx)];
            end
        end
        fields_3d = {'c_s_n_pade', 'c_s_p_pade'};
        for k = 1:numel(fields_3d)
            fn = fields_3d{k};
            if isfield(sim, fn) && isfield(seg, fn)
                sim.(fn) = cat(3, sim.(fn), seg.(fn)(:, :, idx));
            end
        end
    end
    if isfield(sim, 'spme_outputs') && isfield(seg, 'spme_outputs')
        names = fieldnames(sim.spme_outputs);
        for k = 1:numel(names)
            fn = names{k};
            if isfield(seg.spme_outputs, fn)
                sim.spme_outputs.(fn) = [sim.spme_outputs.(fn); seg.spme_outputs.(fn)(idx)];
            end
        end
    end
end

function sim = prune_heavy_history_fields(sim)
    heavy_fields = {'x', 'y', 'z', 'c_s_n', 'c_s_p', 'c_e', 'temperature', ...
                    'c_s_n_pade', 'c_s_p_pade', 'delta_sei', 'delta_li', ...
                    'epsilon_e_n', 'c_ec_sei'};
    for k = 1:numel(heavy_fields)
        fn = heavy_fields{k};
        if isfield(sim, fn)
            sim = rmfield(sim, fn);
        end
    end
end

function sim = apply_voltage_window(sim, vmin_custom, vmax_custom)
    if isempty(vmin_custom) && isempty(vmax_custom)
        return;
    end

    v = sim.voltage(:);
    cut_idx = [];
    if ~isempty(vmin_custom)
        kmin = find(v <= vmin_custom, 1, 'first');
        if ~isempty(kmin)
            cut_idx = kmin;
        end
    end
    if ~isempty(vmax_custom)
        kmax = find(v >= vmax_custom, 1, 'first');
        if isempty(cut_idx)
            cut_idx = kmax;
        elseif ~isempty(kmax)
            cut_idx = min(cut_idx, kmax);
        end
    end
    if isempty(cut_idx) || cut_idx >= numel(v)
        return;
    end

    keep = 1:cut_idx;
    sim.time = sim.time(keep);
    sim.current = sim.current(keep);
    sim.voltage = sim.voltage(keep);
    sim.soc = sim.soc(keep);
    if isfield(sim, 'profile_time'), sim.profile_time = sim.profile_time(keep); end
    if isfield(sim, 'x'), sim.x = sim.x(:, keep); end
    if isfield(sim, 'y'), sim.y = sim.y(:, keep); end
    if isfield(sim, 'z'), sim.z = sim.z(:, keep); end
    if isfield(sim, 'c_s_n'), sim.c_s_n = sim.c_s_n(:, keep); end
    if isfield(sim, 'c_s_p'), sim.c_s_p = sim.c_s_p(:, keep); end
    if isfield(sim, 'c_e'), sim.c_e = sim.c_e(:, keep); end
    if isfield(sim, 'temperature'), sim.temperature = sim.temperature(:, keep); end
    if isfield(sim, 'c_s_n_pade'), sim.c_s_n_pade = sim.c_s_n_pade(:, :, keep); end
    if isfield(sim, 'c_s_p_pade'), sim.c_s_p_pade = sim.c_s_p_pade(:, :, keep); end
    if isfield(sim, 'delta_sei'), sim.delta_sei = sim.delta_sei(keep); end
    if isfield(sim, 'delta_li'), sim.delta_li = sim.delta_li(keep); end
    if isfield(sim, 'epsilon_e_n'), sim.epsilon_e_n = sim.epsilon_e_n(keep); end

    sim.t_last = numel(keep);
    sim.terminated_early = true;
    sim.termination_reason = sprintf('Truncated by voltage window at V=%.6f V.', sim.voltage(end));
end

function update_waitbar_safe(h_wait, idx, total, seg)
    if isempty(h_wait) || ~ishandle(h_wait)
        return;
    end
    frac = min(max(idx / max(total, 1), 0), 1);
    if strcmp(seg.name, 'CC')
        seg_desc = sprintf('CC %.3gC', seg.crate);
    elseif strcmp(seg.name, 'REST')
        seg_desc = sprintf('REST %.3gh', seg.rest_hours);
    else
        seg_desc = seg.name;
    end
    msg = sprintf('Running segments: %d/%d (cycle %d, segment %d: %s)', ...
        idx, total, seg.cycle, seg.segment_in_cycle, seg_desc);
    try
        waitbar(frac, h_wait, msg);
    catch
        % no-op
    end
end

function close_waitbar_safe(h_wait)
    if isempty(h_wait) || ~ishandle(h_wait)
        return;
    end
    try
        close(h_wait);
    catch
        % no-op
    end
end

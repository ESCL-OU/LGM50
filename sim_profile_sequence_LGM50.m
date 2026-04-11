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
    addParameter(parser, 'VEnd', []);
    addParameter(parser, 'Repeat', 1);
    addParameter(parser, 'VMin', []);
    addParameter(parser, 'VMax', []);
    addParameter(parser, 'ShowProgress', true);
    addParameter(parser, 'SPMeOutputs', {});
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

    segments = parse_profile_sequence(profile_sequence);
    segments = repeat_segments(segments, opt.Repeat);
    model_name = upper(char(model_name));
    if ~ismember(model_name, {'DFN', 'SPME', 'ECM'})
        error('model_name must be ''DFN'', ''SPMe'', or ''ECM''.');
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

    seg0 = run_segment(model_name, init_soc, segments(1), opt, struct());
    seg0.sim = apply_voltage_window(seg0.sim, opt.VMin, opt.VMax);
    update_waitbar_safe(h_wait, 1, n_segments_total, segments(1));
    sim = seg0.sim;
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

    for s = 2:numel(segments)
        init_state = terminal_state_from_sim(seg0.sim, model_name);
        init_soc_seg = seg0.sim.soc(end);
        seg0 = run_segment(model_name, init_soc_seg, segments(s), opt, init_state);
        seg0.sim = apply_voltage_window(seg0.sim, opt.VMin, opt.VMax);
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

        sim = append_optional_fields(sim, seg, idx);
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

function out = run_segment(model_name, init_soc, seg, opt, init_state)
    if strcmp(seg.name, 'CC')
        prof = {'CC', seg.crate};
    elseif strcmp(seg.name, 'REST')
        prof = {'REST', seg.rest_hours};
    else
        prof = seg.name;
    end
    tmp_sim = [];
    evalc(['tmp_sim = sim_profile_LGM50(model_name, init_soc, prof, ' ...
        '''Nxn'', opt.Nxn, ''Nxs'', opt.Nxs, ''Nxp'', opt.Nxp, ''Nr'', opt.Nr, ' ...
        '''DeltaT'', opt.DeltaT, ''Tamb'', opt.Tamb, ''ECMOrder'', opt.ECMOrder, ' ...
        '''Crate'', opt.Crate, ''DegradationModel'', opt.DegradationModel, ' ...
        '''ElectrolyteModel'', opt.ElectrolyteModel, ' ...
        '''VEnd'', opt.VEnd, ''InitialState'', init_state, ' ...
        '''SuppressSummary'', true, ''SPMeOutputs'', opt.SPMeOutputs, ''SaveOutput'', false);']);
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
        case 'DFN'
            init_state.x = sim.x(:, end);
            init_state.z = sim.z(:, end);
        case 'ECM'
            init_state.x = sim.x(:, end);
    end
end

function sim = append_optional_fields(sim, seg, idx)
    fields_2d = {'x', 'y', 'z', 'c_s_n', 'c_s_p', 'c_e', 'temperature'};
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

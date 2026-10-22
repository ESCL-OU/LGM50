function out = test_dfn_degradation_spatial_checkpoints(resume_file)
%TEST_DFN_DEGRADATION_SPATIAL_CHECKPOINTS
% Run DFN Safari degradation with Li plating + anode porosity loss and save
% compact spatial checkpoints for later visualization.
%
% This is intended for a long super-computer run. It does NOT make figures.
% Use plot_dfn_degradation_spatial_from_data.m later to visualize the .mat.
%
% Usage:
%   out = test_dfn_degradation_spatial_checkpoints();
%   out = test_dfn_degradation_spatial_checkpoints('results/DFN/dfn_spatial_checkpoint_latest.mat');
%
% Notes:
%   - DFNLplCouplingMode is fixed to 'spme_compatible'.
%   - Saves compact snapshots across x from negative current collector to
%     positive current collector. Variables that only exist in anode are
%     stored as NaN outside the anode for plotting on the same x-axis.

    %% User settings
    degradation_model = 'safari';
    enable_lithium_plating = true;
    enable_anode_porosity_loss = true;
    dfn_lpl_coupling_mode = 'spme_compatible';
    dfn_jacobian_mode = 'padded_analytic';

    init_soc = 0.1;
    sequence1 = {'REST', 0.25, 'CC', -0.2, 'REST', 0.25, 'CC', 0.2};
    sequence2 = {'CC', -1.0, 'CC', 1.0};
    sequence1_repeat = 1;
    sequence2_repeat = 99;
    repeat_count = 20;

    % DFN is much slower than SPMe. Start coarse; increase later if needed.
    delta_t = 1;
    t_amb = 25;
    Nr = 10;
    Nxn = 10;
    Nxs = 5;
    Nxp = 10;

    % Save a compact spatial checkpoint every this many capacity cycles.
    % Set to 1 if you want every cycle, but file size and runtime increase.
    checkpoint_cycle_interval = 10;
    save_each_repeat = true;
    save_latest_checkpoint = true;

    %% Paths
    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    if isempty(base_dir) || strcmp(base_dir, filesep)
        % If this file lives at repo root rather than tests/, use current folder.
        base_dir = pwd;
    end
    addpath(base_dir);
    addpath(fullfile(base_dir, 'DFN'));
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    results_dir = fullfile(base_dir, 'results', 'DFN');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end
    run_ts = datestr(now, 'yyyymmdd_HHMMSS');
    output_file = fullfile(results_dir, sprintf('dfn_degradation_spatial_checkpoints_%s.mat', run_ts));
    latest_file = fullfile(results_dir, 'dfn_degradation_spatial_checkpoint_latest.mat');

    %% Build sequence
    block_sequence = [repmat(sequence1, 1, sequence1_repeat), ...
                      repmat(sequence2, 1, sequence2_repeat)];
    segments_base = parse_profile_sequence_local(block_sequence);
    segments = repeat_segments_local(segments_base, repeat_count);
    capacity_cycle_by_segment = assign_capacity_cycle_index_local(segments);

    %% Resume support
    start_segment = 1;
    initial_state = struct();
    global_time_offset = 0;
    charged_ah_cycle = 0;
    discharged_ah_cycle = 0;
    current_capacity_cycle = capacity_cycle_by_segment(1);
    snapshots = struct([]);
    segment_log = struct([]);
    completed_segments = 0;
    init_soc_run = init_soc;
    resumed_from_checkpoint = false;

    if nargin >= 1 && ~isempty(resume_file)
        loaded = load(resume_file, 'checkpoint');
        if ~isfield(loaded, 'checkpoint')
            error('Resume file does not contain variable checkpoint: %s', resume_file);
        end
        checkpoint = loaded.checkpoint;
        if ~strcmpi(checkpoint.model, 'DFN')
            error('Resume file is not a DFN checkpoint.');
        end
        snapshots = checkpoint.snapshots;
        segment_log = checkpoint.segment_log;
        completed_segments = checkpoint.completed_segments;
        start_segment = completed_segments + 1;
        initial_state = checkpoint.terminal_state;
        global_time_offset = checkpoint.global_time_offset;
        charged_ah_cycle = checkpoint.charged_ah_cycle;
        discharged_ah_cycle = checkpoint.discharged_ah_cycle;
        current_capacity_cycle = checkpoint.current_capacity_cycle;
        init_soc_run = checkpoint.terminal_soc;
        resumed_from_checkpoint = true;
        fprintf('Resuming from %s at segment %d / %d.\n', resume_file, start_segment, numel(segments));
    end

    settings = struct();
    settings.degradation_model = degradation_model;
    settings.enable_lithium_plating = enable_lithium_plating;
    settings.enable_anode_porosity_loss = enable_anode_porosity_loss;
    settings.dfn_lpl_coupling_mode = dfn_lpl_coupling_mode;
    settings.dfn_jacobian_mode = dfn_jacobian_mode;
    settings.init_soc = init_soc;
    settings.sequence1 = {sequence1};
    settings.sequence2 = {sequence2};
    settings.sequence1_repeat = sequence1_repeat;
    settings.sequence2_repeat = sequence2_repeat;
    settings.repeat_count = repeat_count;
    settings.delta_t = delta_t;
    settings.t_amb = t_amb;
    settings.Nr = Nr;
    settings.Nxn = Nxn;
    settings.Nxs = Nxs;
    settings.Nxp = Nxp;
    settings.checkpoint_cycle_interval = checkpoint_cycle_interval;
    settings.save_each_repeat = save_each_repeat;
    settings.run_timestamp = run_ts;

    fprintf('DFN spatial checkpoint run started.\n');
    fprintf('Output file: %s\n', output_file);
    fprintf('Mode: %s\n', dfn_lpl_coupling_mode);
    fprintf('Jacobian mode: %s\n', dfn_jacobian_mode);
    fprintf('Total segments: %d\n', numel(segments));

    %% Main segment loop
    memory_light = true;
    show_progress = true;
    enable_diagnostics = false;
    diagnostics_file = '';
    diagnostic_context = struct('test_name', mfilename, 'run_timestamp', run_ts);
    checkpoint_file = latest_file;
    callback = @(seg,seg_idx,seg_info,cycle_idx) ...
        dfn_segment_checkpoint_callback(seg, seg_idx, seg_info, cycle_idx, ...
        output_file, latest_file, settings, checkpoint_cycle_interval, sequence1_repeat, sequence2_repeat);
    
    sim = sim_profile_sequence_LGM50('DFN', init_soc_run, block_sequence, ...
        'Repeat', repeat_count, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', degradation_model, ...
        'EnableLithiumPlating', enable_lithium_plating, ...
        'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
        'DFNLplCouplingMode', dfn_lpl_coupling_mode, ...
        'DFNJacobianMode', dfn_jacobian_mode, ...
        'SaveOutput', false, ...
        'ShowProgress', show_progress, ...
        'MemoryLight', memory_light, ...
        'InitialState', initial_state, ...
        'CheckpointFile', checkpoint_file, ...
        'EnableDiagnostics', enable_diagnostics, ...
        'DiagnosticsFile', diagnostics_file, ...
        'DiagnosticContext', diagnostic_context, ...
        'SegmentCallback', callback);

    loaded = load(latest_file, 'checkpoint');
    out = loaded.checkpoint;
    fprintf('DFN spatial checkpoint run completed.\nSaved: %s\n', output_file);
end

function dfn_segment_checkpoint_callback(seg, seg_idx, seg_info, cycle_idx, ...
        output_file, latest_file, settings, checkpoint_cycle_interval, ...
        sequence1_repeat, sequence2_repeat)

    persistent snapshots segment_log

    if isempty(snapshots)
        snapshots = struct([]);
    end
    if isempty(segment_log)
        segment_log = struct([]);
    end

    charged_ah = trapz(seg.time(:), max(-seg.current(:),0)) / 3600;
    discharged_ah = trapz(seg.time(:), max(seg.current(:),0)) / 3600;

    log_entry = struct();
    log_entry.segment = seg_idx;
    log_entry.capacity_cycle = cycle_idx;
    log_entry.name = seg_info.name;
    log_entry.crate = seg_info.crate;
    log_entry.rest_hours = seg_info.rest_hours;
    log_entry.terminal_soc = seg.soc(end);
    log_entry.terminal_voltage = seg.voltage(end);
    log_entry.charged_Ah = charged_ah;
    log_entry.discharged_Ah = discharged_ah;
    log_entry.net_Ah = discharged_ah - charged_ah;
    log_entry.terminated_early = seg.terminated_early;
    log_entry.termination_reason = seg.termination_reason;

    if isempty(segment_log)
        segment_log = log_entry;
    else
        segment_log(end+1) = log_entry; %#ok<AGROW>
    end

    [snapshot_tag, should_snapshot] = classify_dfn_checkpoint(seg_info, cycle_idx, ...
        sequence1_repeat, sequence2_repeat);

    if should_snapshot
        snap = make_dfn_spatial_snapshot_local(seg, cycle_idx, seg_idx, ...
            seg.time(end), charged_ah, discharged_ah);
        snap.tag = snapshot_tag;
        snap.cycle_idx = cycle_idx;
        snap.segment_name = seg_info.name;
        snap.segment_crate = seg_info.crate;

        if isempty(snapshots)
            snapshots = snap;
        else
            snapshots(end+1) = snap; %#ok<AGROW>
        end
    end

    checkpoint = struct();
    checkpoint.model = 'DFN';
    checkpoint.saved_at = datestr(now, 31);
    checkpoint.settings = settings;
    checkpoint.snapshots = snapshots;
    checkpoint.segment_log = segment_log;
    checkpoint.completed_segments = seg_idx;
    checkpoint.terminal_state = terminal_state_from_sim_local(seg);
    checkpoint.terminal_soc = seg.soc(end);
    checkpoint.terminal_voltage = seg.voltage(end);
    checkpoint.latest_file = latest_file;
    checkpoint.output_file = output_file;

    save(latest_file, 'checkpoint', '-v7.3');

    if should_snapshot
        save(output_file, 'checkpoint', '-v7.3');
        fprintf('  Saved DFN spatial snapshot at segment %d, cycle %d.\n', seg_idx, cycle_idx);
        drawnow;
    end
end

function snapshot = make_dfn_spatial_snapshot_local(sim, cap_cycle, segment_idx, global_time_s, charged_ah, discharged_ah)
    p = sim.p;
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    Nce = p.Nx - 3;

    d = sim.dfn_outputs;
    k = numel(sim.time);
    x_info = dfn_x_axis_local(p);

    c_e = nan(Nce, 1);
    if isfield(sim, 'c_e') && ~isempty(sim.c_e)
        c_e = sim.c_e(:, end);
    end

    LDP_all = nan(Nce, 1);
    LDP_total_all = nan(Nce, 1);
    j_lpl_all = nan(Nce, 1);
    j_sei_all = nan(Nce, 1);
    delta_sei_all = nan(Nce, 1);
    delta_li_all = nan(Nce, 1);
    epsilon_all = nan(Nce, 1);

    LDP_all(1:Nn) = d.LDP_n(:, end);
    if isfield(d, 'LDP_total_n')
        LDP_total_all(1:Nn) = d.LDP_total_n(:, end);
    end
    j_lpl_all(1:Nn) = d.j_lpl_n(:, end);
    j_sei_all(1:Nn) = d.j_sei_n(:, end);
    delta_sei_all(1:Nn) = d.delta_sei_n(:, end);
    delta_li_all(1:Nn) = d.delta_li_n(:, end);
    epsilon_all(1:Nn) = d.epsilon_e_n_dyn(:, end);
    epsilon_all(Nn+1:Nn+Ns) = p.epsilon_e_s;
    epsilon_all(Nn+Ns+1:end) = p.epsilon_e_p;

    snapshot = struct();
    snapshot.capacity_cycle = cap_cycle;
    snapshot.segment_idx = segment_idx;
    snapshot.global_time_s = global_time_s;
    snapshot.local_time_s = sim.time(end);
    snapshot.current_A = sim.current(end);
    snapshot.voltage_V = sim.voltage(end);
    snapshot.soc = sim.soc(end);
    snapshot.charged_Ah_in_cycle = charged_ah;
    snapshot.discharged_Ah_in_cycle = discharged_ah;
    snapshot.net_Ah_in_cycle = discharged_ah - charged_ah;
    snapshot.x_norm = x_info.x_norm;
    snapshot.x_region = x_info.region_id;
    snapshot.region_names = x_info.region_names;
    snapshot.x_boundaries_norm = x_info.boundaries_norm;
    snapshot.c_e = c_e;
    snapshot.LDP = LDP_all;
    snapshot.LDP_total = LDP_total_all;
    snapshot.j_lpl = j_lpl_all;
    snapshot.j_sei = j_sei_all;
    snapshot.delta_sei = delta_sei_all;
    snapshot.delta_li = delta_li_all;
    snapshot.epsilon_e = epsilon_all;
    snapshot.separator_side = struct( ...
        'LDP', d.LDP_sep_n(end), ...
        'j_lpl', d.j_lpl_sep_n(end), ...
        'delta_li', d.delta_li_sep_n(end), ...
        'epsilon_e', d.epsilon_e_sep_n(end));
    if isfield(d, 'LDP_total_sep_n')
        snapshot.separator_side.LDP_total = d.LDP_total_sep_n(end);
    end
    snapshot.discretization = sim.discretization;
end

function x_info = dfn_x_axis_local(p)
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    Ln = p.L_n;
    Ls = p.L_s;
    Lp = p.L_p;
    Ltot = Ln + Ls + Lp;

    x_n = linspace(0, Ln, Nn+2).';
    x_n = x_n(2:end-1);
    x_s = linspace(Ln, Ln+Ls, Ns+2).';
    x_s = x_s(2:end-1);
    x_p = linspace(Ln+Ls, Ltot, Np+2).';
    x_p = x_p(2:end-1);
    x = [x_n; x_s; x_p] / Ltot;

    x_info = struct();
    x_info.x_norm = x;
    x_info.region_id = [ones(Nn,1); 2*ones(Ns,1); 3*ones(Np,1)];
    x_info.region_names = {'Anode', 'Separator', 'Cathode'};
    x_info.boundaries_norm = [0, Ln/Ltot, (Ln+Ls)/Ltot, 1];
end

function checkpoint = make_checkpoint_local(settings, snapshots, segment_log, completed_segments, total_segments, terminal_state, terminal_soc, global_time_offset, charged_ah_cycle, discharged_ah_cycle, current_capacity_cycle, output_file, latest_file, resumed_from_checkpoint)
    checkpoint = struct();
    checkpoint.model = 'DFN';
    checkpoint.saved_at = datestr(now, 31);
    checkpoint.settings = settings;
    checkpoint.snapshots = snapshots;
    checkpoint.segment_log = segment_log;
    checkpoint.completed_segments = completed_segments;
    checkpoint.total_segments = total_segments;
    checkpoint.completed_fraction = completed_segments / max(total_segments, 1);
    checkpoint.terminal_state = terminal_state;
    checkpoint.terminal_soc = terminal_soc;
    checkpoint.global_time_offset = global_time_offset;
    checkpoint.charged_ah_cycle = charged_ah_cycle;
    checkpoint.discharged_ah_cycle = discharged_ah_cycle;
    checkpoint.current_capacity_cycle = current_capacity_cycle;
    checkpoint.output_file = output_file;
    checkpoint.latest_file = latest_file;
    checkpoint.resumed_from_checkpoint = resumed_from_checkpoint;
end

function segments = parse_profile_sequence_local(profile_sequence)
    if ischar(profile_sequence) || isstring(profile_sequence)
        segments = struct('name', upper(char(profile_sequence)), 'crate', NaN, 'rest_hours', NaN, 'cycle', 1, 'segment_in_cycle', 1);
        return;
    end
    if ~iscell(profile_sequence)
        error('profile_sequence must be char/string or cell.');
    end
    segments = struct('name', {}, 'crate', {}, 'rest_hours', {}, 'cycle', {}, 'segment_in_cycle', {});
    i = 1;
    seg_in_cycle = 1;
    cycle = 1;
    last_nonzero_dir = 0;
    while i <= numel(profile_sequence)
        name = upper(char(profile_sequence{i}));
        if strcmp(name, 'CC')
            val = profile_sequence{i+1};
            dir = sign(val);
            if last_nonzero_dir ~= 0 && dir ~= 0 && dir ~= last_nonzero_dir
                cycle = cycle + 0.5;
            end
            if dir ~= 0
                last_nonzero_dir = dir;
            end
            segments(end+1) = struct('name', 'CC', 'crate', val, 'rest_hours', NaN, ...
                'cycle', ceil(cycle/2), 'segment_in_cycle', seg_in_cycle); %#ok<AGROW>
            i = i + 2;
        elseif strcmp(name, 'REST')
            val = profile_sequence{i+1};
            segments(end+1) = struct('name', 'REST', 'crate', 0, 'rest_hours', val, ...
                'cycle', ceil(max(cycle,1)/2), 'segment_in_cycle', seg_in_cycle); %#ok<AGROW>
            i = i + 2;
        else
            segments(end+1) = struct('name', name, 'crate', NaN, 'rest_hours', NaN, ...
                'cycle', ceil(max(cycle,1)/2), 'segment_in_cycle', seg_in_cycle); %#ok<AGROW>
            i = i + 1;
        end
        seg_in_cycle = seg_in_cycle + 1;
    end
end

function segments_rep = repeat_segments_local(segments_base, repeat_count)
    segments_rep = struct('name', {}, 'crate', {}, 'rest_hours', {}, 'cycle', {}, 'segment_in_cycle', {});
    n = numel(segments_base);
    for r = 1:repeat_count
        for s = 1:n
            segments_rep(end+1) = segments_base(s); %#ok<AGROW>
            segments_rep(end).cycle = (r-1) * max([segments_base.cycle]) + segments_base(s).cycle;
        end
    end
end

function cap_cycle = assign_capacity_cycle_index_local(segments)
    cap_cycle = [segments.cycle].';
end

function prof = segment_to_profile_def_local(seg)
    if strcmp(seg.name, 'CC')
        prof = {'CC', seg.crate};
    elseif strcmp(seg.name, 'REST')
        prof = {'REST', seg.rest_hours};
    else
        prof = seg.name;
    end
end

function label = segment_label_local(seg)
    if strcmp(seg.name, 'CC')
        label = sprintf('CC %.3gC', seg.crate);
    elseif strcmp(seg.name, 'REST')
        label = sprintf('REST %.3gh', seg.rest_hours);
    else
        label = seg.name;
    end
end

function [charged_ah, discharged_ah] = compute_charge_discharge_ah_local(time_s, current_a)
    if isempty(time_s) || numel(time_s) < 2
        charged_ah = 0;
        discharged_ah = 0;
        return;
    end
    t = time_s(:);
    i = current_a(:);
    n = min(numel(t), numel(i));
    t = t(1:n);
    i = i(1:n);
    charged_ah = trapz(t, max(-i,0)) / 3600;
    discharged_ah = trapz(t, max(i,0)) / 3600;
end

function init_state = terminal_state_from_sim_local(sim)
    init_state = struct();
    init_state.x = sim.x(:, end);
    init_state.z = sim.z(:, end);
    init_state.p = sim.p;
end

function [tag, tf] = classify_dfn_checkpoint(seg_info, cycle_idx, sequence1_repeat, sequence2_repeat)

    tag = '';
    tf = false;

    if ~strcmp(seg_info.name, 'CC')
        return;
    end

    crate = seg_info.crate;

    if abs(crate + 0.2) < 1e-12
        tag = '0p2C_charge';
        tf = true;
        return;
    elseif abs(crate - 0.2) < 1e-12
        tag = '0p2C_discharge';
        tf = true;
        return;
    end

    if abs(crate + 1.0) < 1e-12
        base_tag = '1C_charge';
    elseif abs(crate - 1.0) < 1e-12
        base_tag = '1C_discharge';
    else
        return;
    end

    first_cycle = sequence1_repeat + 1;
    middle_cycle = sequence1_repeat + ceil(sequence2_repeat/2);
    last_cycle = sequence1_repeat + sequence2_repeat;

    if ismember(cycle_idx, [first_cycle, middle_cycle, last_cycle])
        tag = base_tag;
        tf = true;
    end
end

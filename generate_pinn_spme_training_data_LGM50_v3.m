function generate_pinn_spme_training_data_LGM50_v3(checkpoint_file)
%GENERATE_PINN_SPME_TRAINING_DATA_LGM50
% Generate PINN-SPMe training/validation datasets using the LGM50 SPMe code.
% Optional input:
%   checkpoint_file : .mat checkpoint from test_spme_degradation_addons_sequence_plot.
%                     When provided, only degradation states/parameters
%                     are copied; each dataset still uses the job init_soc.
%
% This replaces the older SPM_FDM hand-tuned current-profile scripts.
% It uses sim_profile_LGM50 / sim_profile_sequence_LGM50 and lets the
% simulator terminate dynamic profiles at voltage limits instead of manually
% trimming/scaling current profiles.
%
% Output:
%   results/PINN_SPMe_dataset/*.mat
%   results/PINN_SPMe_dataset/*.csv
%
% Current sign convention in LGM50:
%   + current = discharge
%   - current = charge
%
% Recommended first use:
%   1) Put this file in the root of the LGM50 repo or in the same folder
%      that contains sim_profile_LGM50.m and sim_profile_sequence_LGM50.m.
%   2) Run:
%          generate_pinn_spme_training_data_LGM50
%
% Notes for PINN-SPMe:
%   MATLAB SPMe state order is [c_s_n; c_s_p; c_e].
%   PINN-SPMe usually uses [c_s_n; c_e; c_s_p].
%   This script saves both:
%       sim.x_matlab_order
%       data.x_pinn_order

if nargin < 1
    checkpoint_file = '';
end
checkpoint_file = text_scalar(checkpoint_file);
clc;
fprintf('generate_pinn_spme_training_data_LGM50_v3 version: path-safe-degraded-params-fresh-soc-20260702c\n');


%% ---------------- User settings ----------------
opt = struct();

% Model/discretization should match PINN-SPMe config.
opt.model_name = 'SPMe';
opt.Nr  = 15;
opt.Nxn = 10;
opt.Nxs = 5;
opt.Nxp = 10;
opt.DeltaT = 1.0;          % [s], target PINN sampling time
opt.Tamb = 25;             % [degC]
opt.ElectrolyteModel = 'fvm';
opt.DegradationModel = 'none';
opt.UseElectrolyteActivity = false;

% Optional degraded-cell generation from checkpoint.
% Mode 'degraded_params_fresh_soc' means: keep the same init_soc logic as
% this generator, but copy degradation states such as delta_sei, delta_li,
% epsilon_e_n, and Safari c_ec_sei from the checkpoint.
opt.CheckpointFile = text_scalar(checkpoint_file);
opt.CheckpointMode = 'degraded_params_fresh_soc';
opt.checkpoint_info = struct();

% High-to-low SoC dataset.
opt.init_soc_discharge = 0.90;
opt.VMin = 2.50;           % stop discharge/dynamic when V <= VMin
opt.VMax = 4.20;           % safety upper limit

% Optional low-to-high charge dataset if needed later.
opt.init_soc_charge = 0.20;
opt.charge_VEnd = 4.20;

% Output.
opt.output_dir = text_scalar(fullfile(pwd, 'results', 'PINN_SPMe_dataset'));
opt.write_csv = true;
opt.write_mat = true;

% Profiles for training/validation.
% Dynamic profiles will be repeated internally until VEnd is reached.
% Use BJDST/FUDS/US06 for dynamic training, DST for validation.
jobs = [
    make_job("train_BJDST", "dynamic", "BJDST", opt.init_soc_discharge, opt.VMin)
    make_job("train_FUDS",  "dynamic", "FUDS",  opt.init_soc_discharge, opt.VMin)
    make_job("train_US06",  "dynamic", "US06",  opt.init_soc_discharge, opt.VMin)
    make_job("val_DST",     "dynamic", "DST",   opt.init_soc_discharge, opt.VMin)
    make_job("debug_CC_1C_discharge", "cc", "CC", opt.init_soc_discharge, opt.VMin, 1.0)
];

add_lgm50_paths();

% Prepare optional checkpoint-derived degraded initial condition.
% This intentionally does NOT use checkpoint.terminal_soc or terminal_state.x
% directly for solid/electrolyte states, because the PINN datasets should
% start from the same requested init_soc as fresh-cell generation.
if ~isempty(opt.CheckpointFile)
    [opt, jobs] = configure_checkpoint_degraded_cell(opt, jobs);
end
opt.output_dir = text_scalar(opt.output_dir);
if ~isfolder(opt.output_dir)
    mkdir(opt.output_dir);
end

%% ---------------- Run jobs ----------------
summary = struct([]);
for j = 1:numel(jobs)
    job = jobs(j);
    fprintf('\n============================================================\n');
    fprintf('Generating %s (%s, profile=%s, init_soc=%.3f)\n', ...
        text_scalar(job.name), text_scalar(job.kind), text_scalar(job.profile), job.init_soc);
    fprintf('============================================================\n');

    sim = run_one_dataset(job, opt);
    data = convert_sim_to_pinn_dataset(sim, opt);

    out_base = text_scalar(fullfile(text_scalar(opt.output_dir), text_scalar(job.name)));
    mat_file = sprintf('%s.mat', out_base);
    csv_file = sprintf('%s.csv', out_base);
    if opt.write_mat
        save(mat_file, 'sim', 'data', 'opt', 'job', '-v7.3');
    end
    if opt.write_csv
        write_core_csv(csv_file, data);
    end

    summary(j).name = text_scalar(job.name); %#ok<AGROW>
    summary(j).n = numel(data.time);
    summary(j).duration_s = data.time(end) - data.time(1);
    summary(j).init_soc = job.init_soc;
    summary(j).final_voltage = data.voltage(end);
    summary(j).final_soc = data.soc(end);
    summary(j).min_voltage = min(data.voltage);
    summary(j).max_voltage = max(data.voltage);
    summary(j).terminated_early = sim.terminated_early;
    summary(j).termination_reason = sim.termination_reason;
    fprintf('Saved: %s.mat / %s.csv\n', text_scalar(out_base), text_scalar(out_base));
    fprintf('Duration %.1f s, N=%d, V %.3f -> %.3f, SoC %.3f -> %.3f\n', ...
        summary(j).duration_s, summary(j).n, data.voltage(1), data.voltage(end), ...
        data.soc(1), data.soc(end));
end

summary_mat_file = text_scalar(fullfile(opt.output_dir, 'dataset_summary.mat'));
summary_csv_file = text_scalar(fullfile(opt.output_dir, 'dataset_summary.csv'));
save(summary_mat_file, 'summary', 'opt', 'jobs');
write_summary_csv(summary_csv_file, summary);

fprintf('\nDone. Output folder:\n  %s\n', text_scalar(opt.output_dir));

end

%% ========================================================================
function job = make_job(name, kind, profile, init_soc, v_end, crate)
    if nargin < 6
        crate = NaN;
    end
    job = struct();
    job.name = text_scalar(name);
    job.kind = text_scalar(kind);
    job.profile = text_scalar(profile);
    job.init_soc = init_soc;
    job.v_end = v_end;
    job.crate = crate;
end

function add_lgm50_paths()
    this_dir = fileparts(mfilename('fullpath'));
    candidate_dirs = {this_dir, fullfile(this_dir, 'SPMe'), fullfile(this_dir, 'param'), ...
                      fullfile(this_dir, 'DFN'), fullfile(this_dir, 'ECM'), fullfile(this_dir, 'ECM 2')};
    for ii = 1:numel(candidate_dirs)
        if isfolder(candidate_dirs{ii})
            addpath(candidate_dirs{ii});
        end
    end
    % Do not call startup.m here; some startup files add optional missing folders
    % such as ECM and create noisy warnings. Candidate folders above are enough.

end


function [opt, jobs] = configure_checkpoint_degraded_cell(opt, jobs)
    opt.CheckpointFile = text_scalar(opt.CheckpointFile);
    if ~isfile(opt.CheckpointFile)
        error('CheckpointFile does not exist: %s', opt.CheckpointFile);
    end
    S = load(opt.CheckpointFile, 'checkpoint');
    if ~isfield(S, 'checkpoint')
        error('Checkpoint file does not contain variable "checkpoint": %s', opt.CheckpointFile);
    end
    ckpt = S.checkpoint;
    p_ckpt = get_checkpoint_params(ckpt);

    if isfield(p_ckpt, 'degradation_model') && ~isempty(p_ckpt.degradation_model)
        opt.DegradationModel = lower(char(p_ckpt.degradation_model));
    end
    if isfield(p_ckpt, 'electrolyte_model') && ~isempty(p_ckpt.electrolyte_model)
        opt.ElectrolyteModel = lower(char(p_ckpt.electrolyte_model));
    end
    if isfield(p_ckpt, 'use_electrolyte_activity') && ~isempty(p_ckpt.use_electrolyte_activity)
        opt.UseElectrolyteActivity = logical(p_ckpt.use_electrolyte_activity);
    end

    deg = extract_degradation_state_from_checkpoint(ckpt, p_ckpt);
    opt.checkpoint_info.file = text_scalar(opt.CheckpointFile);
    opt.checkpoint_info.mode = opt.CheckpointMode;
    opt.checkpoint_info.terminal_soc_ignored_for_init = getfield_or(ckpt, 'terminal_soc', NaN);
    opt.checkpoint_info.degradation_model = opt.DegradationModel;
    opt.checkpoint_info.electrolyte_model = opt.ElectrolyteModel;
    opt.checkpoint_info.degraded_state = deg;

    % Keep static compatibility parameters from checkpoint. Do not override
    % csn0/csp0/x0, because those must be rebuilt from each job init_soc.
    opt.SPMeParamOverrides = struct();
    copy_names = {'Nsei','epsilon_sei','delta_sei_min', ...
                  'enable_lithium_plating','enable_anode_porosity_loss', ...
                  'i0_lpl','alph_lpl','M_li','rho_li', ...
                  'kappa_sei','M_sei','rho_sei','w_sei', ...
                  'D_ec','c0_ec','k_sei','beta_sei'};
    for i = 1:numel(copy_names)
        nm = copy_names{i};
        if isfield(p_ckpt, nm)
            opt.SPMeParamOverrides.(nm) = p_ckpt.(nm);
        end
    end

    % Add checkpoint tag to output folder so fresh and degraded datasets do not collide.
    [~, ckpt_base, ~] = fileparts(text_scalar(opt.CheckpointFile));
    safe_base = regexprep(text_scalar(ckpt_base), '[^\w\-]', '_');
    opt.output_dir = text_scalar(fullfile(text_scalar(opt.output_dir), ['degraded_from_' text_scalar(safe_base)]));

    for j = 1:numel(jobs)
        jobs(j).checkpoint_file = text_scalar(opt.CheckpointFile); %#ok<AGROW>
        jobs(j).checkpoint_mode = opt.CheckpointMode; %#ok<AGROW>
    end

    fprintf('\nUsing degraded parameters/states from checkpoint:\n  %s\n', opt.CheckpointFile);
    fprintf('  checkpoint terminal_soc = %.6g is ignored for dataset initialization.\n', opt.checkpoint_info.terminal_soc_ignored_for_init);
    fprintf('  Each job still starts from its configured init_soc.\n');
    fprintf('  Degradation model = %s, electrolyte model = %s\n', opt.DegradationModel, opt.ElectrolyteModel);
    fprintf('  delta_sei = %.6g m, delta_li = %.6g m, epsilon_e_n = %.6g\n', ...
        deg.delta_sei, deg.delta_li, deg.epsilon_e_n);
end

function p = get_checkpoint_params(ckpt)
    p = struct();
    if isfield(ckpt, 'terminal_state') && isfield(ckpt.terminal_state, 'p') && ~isempty(ckpt.terminal_state.p)
        p = ckpt.terminal_state.p;
    elseif isfield(ckpt, 'last_segment') && isfield(ckpt.last_segment, 'p') && ~isempty(ckpt.last_segment.p)
        p = ckpt.last_segment.p;
    elseif isfield(ckpt, 'sim') && isfield(ckpt.sim, 'p') && ~isempty(ckpt.sim.p)
        p = ckpt.sim.p;
    else
        error('Could not find checkpoint parameters in terminal_state.p, last_segment.p, or sim.p.');
    end
end

function deg = extract_degradation_state_from_checkpoint(ckpt, p)
    if ~isfield(ckpt, 'terminal_state') || ~isfield(ckpt.terminal_state, 'x') || isempty(ckpt.terminal_state.x)
        error('Checkpoint is missing terminal_state.x.');
    end
    x = ckpt.terminal_state.x(:);
    idx = degradation_state_indices(p);

    deg = struct();
    deg.delta_sei = 0;
    deg.delta_li = 0;
    deg.epsilon_e_n = getfield_or(p, 'epsilon_e_n', NaN);
    deg.c_ec_sei = [];

    if ~isnan(idx.delta_sei)
        deg.delta_sei = x(idx.delta_sei);
    end
    if ~isnan(idx.delta_li)
        deg.delta_li = x(idx.delta_li);
    end
    if ~isnan(idx.epsilon_e_n)
        deg.epsilon_e_n = x(idx.epsilon_e_n);
    end
    if ~isempty(idx.c_ec_sei)
        deg.c_ec_sei = x(idx.c_ec_sei);
    end
end

function idx = degradation_state_indices(p)
    Nr = p.Nr;
    Nx = p.Nx;
    base = 2 * (Nr - 1) + (Nx - 3);
    k = base;

    model = 'none';
    if isfield(p, 'degradation_model') && ~isempty(p.degradation_model)
        model = lower(text_scalar(p.degradation_model));
    end
    has_sei = ~strcmpi(model, 'none');
    has_plating = isfield(p, 'enable_lithium_plating') && logical(p.enable_lithium_plating);
    has_porosity = isfield(p, 'enable_anode_porosity_loss') && logical(p.enable_anode_porosity_loss);
    is_safari = strcmpi(model, 'safari');

    idx = struct('delta_sei', NaN, 'delta_li', NaN, 'epsilon_e_n', NaN, 'c_ec_sei', []);
    if has_sei
        k = k + 1;
        idx.delta_sei = k;
    end
    if has_plating
        k = k + 1;
        idx.delta_li = k;
    end
    if has_porosity
        k = k + 1;
        idx.epsilon_e_n = k;
    end
    if is_safari
        nsei = getfield_or(p, 'Nsei', 10);
        idx.c_ec_sei = (k+1):(k+nsei);
    end
end

function initial_state = make_job_initial_state(job, opt)
    initial_state = struct();
    if ~isfield(opt, 'checkpoint_info') || ~isfield(opt.checkpoint_info, 'degraded_state') || isempty(fieldnames(opt.checkpoint_info))
        return;
    end

    % Build a fresh SoC-consistent SPMe x0 for this job, then append/copy
    % checkpoint degradation states. This mirrors the simulator's own fresh
    % initialization logic, but uses the checkpoint degradation condition.
    % Suppress load_SPMe_params console output because it always prints
    % fresh-cell delta_sei0 = 0 before we inject checkpoint degradation states.
    [~, p] = evalc('load_SPMe_params(job.init_soc, opt.Tamb, opt.Nr, opt.DeltaT, opt.Nxn, opt.Nxs, opt.Nxp, true, opt.DegradationModel, opt.ElectrolyteModel, opt.UseElectrolyteActivity, get_override_or_default(opt, ''enable_lithium_plating'', []), get_override_or_default(opt, ''enable_anode_porosity_loss'', []))');
    p = apply_local_param_overrides(p, get_job_param_overrides(opt));

    calc = calculate_equilibrium_voltage_from_soc(job.init_soc);
    [csn0, csp0] = init_cs(p, calc.Veq);
    p.csn0 = csn0;
    p.csp0 = csp0;

    x0 = [csn0 * ones(p.Nr - 1, 1); ...
          csp0 * ones(p.Nr - 1, 1); ...
          p.c_e * ones(p.Nx - 3, 1)];

    deg = opt.checkpoint_info.degraded_state;
    if ~strcmpi(text_scalar(opt.DegradationModel), 'none')
        x0 = [x0; deg.delta_sei];
    end
    if isfield(p, 'enable_lithium_plating') && logical(p.enable_lithium_plating)
        x0 = [x0; deg.delta_li];
    end
    if isfield(p, 'enable_anode_porosity_loss') && logical(p.enable_anode_porosity_loss)
        x0 = [x0; deg.epsilon_e_n];
    end
    if strcmpi(text_scalar(opt.DegradationModel), 'safari')
        nsei = getfield_or(p, 'Nsei', 10);
        if isfield(deg, 'c_ec_sei') && numel(deg.c_ec_sei) == nsei
            c_ec = deg.c_ec_sei(:);
        else
            c_ec = p.epsilon_sei * p.c0_ec * ones(nsei, 1);
        end
        x0 = [x0; c_ec];
    end

    p.x0 = x0;
    initial_state.x = x0;
    initial_state.p = p;
    initial_state.degraded_state_source = opt.checkpoint_info.file;
    initial_state.init_soc_rebuilt_from_job = job.init_soc;
    initial_state.delta_sei0_from_checkpoint = deg.delta_sei;
    initial_state.delta_li0_from_checkpoint = deg.delta_li;
    initial_state.epsilon_e_n0_from_checkpoint = deg.epsilon_e_n;
end

function overrides = get_job_param_overrides(opt)
    overrides = struct();
    if isfield(opt, 'SPMeParamOverrides') && isstruct(opt.SPMeParamOverrides)
        overrides = opt.SPMeParamOverrides;
    end
end

function val = get_override_or_default(opt, name, default_val)
    val = default_val;
    if isfield(opt, 'SPMeParamOverrides') && isfield(opt.SPMeParamOverrides, name)
        val = opt.SPMeParamOverrides.(name);
    end
end

function p = apply_local_param_overrides(p, overrides)
    if isempty(overrides) || ~isstruct(overrides)
        return;
    end
    names = fieldnames(overrides);
    for i = 1:numel(names)
        p.(names{i}) = overrides.(names{i});
    end
end

function txt = text_scalar(x)
    % Return a MATLAB char row vector from char/string/cellstr-like input.
    if iscell(x)
        if isempty(x)
            txt = '';
        else
            txt = text_scalar(x{1});
        end
    elseif isstring(x)
        if isempty(x)
            txt = '';
        else
            txt = char(x(1));
        end
    elseif ischar(x)
        txt = x;
    else
        txt = char(string(x));
    end
end

function val = getfield_or(s, name, default_val)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default_val;
    end
end

function sim = run_one_dataset(job, opt)
    init_state = make_job_initial_state(job, opt);
    if isfield(init_state, 'delta_sei0_from_checkpoint')
        fprintf('Injected checkpoint degradation at fresh init_soc %.3f: delta_sei=%.6g m (%.6f um), delta_li=%.6g m, epsilon_e_n=%.6g\n', ...
            job.init_soc, init_state.delta_sei0_from_checkpoint, init_state.delta_sei0_from_checkpoint*1e6, ...
            init_state.delta_li0_from_checkpoint, init_state.epsilon_e_n0_from_checkpoint);
    end
    common_args = { ...
        'Nxn', opt.Nxn, 'Nxs', opt.Nxs, 'Nxp', opt.Nxp, 'Nr', opt.Nr, ...
        'DeltaT', opt.DeltaT, ...
        'Tamb', opt.Tamb, ...
        'DegradationModel', opt.DegradationModel, ...
        'ElectrolyteModel', opt.ElectrolyteModel, ...
        'UseElectrolyteActivity', opt.UseElectrolyteActivity, ...
        'SPMeParamOverrides', get_job_param_overrides(opt), ...
        'InitialState', init_state, ...
        'SPMeOutputs', {'V_electrolytePolar','ce0n','cens','cesp','ce0p', ...
                        'cen_bar','ces_bar','cep_bar','cs_ave_n','cs_ave_p'}, ...
        'SaveOutput', false, ...
        'SuppressSummary', isfield(init_state, 'delta_sei0_from_checkpoint')};

    switch lower(text_scalar(job.kind))
        case 'dynamic'
            % Dynamic profile is repeated by sim_profile_LGM50 until VEnd.
            sim = sim_profile_LGM50(opt.model_name, job.init_soc, text_scalar(job.profile), ...
                common_args{:}, ...
                'VEnd', job.v_end);

        case 'cc'
            sim = sim_profile_LGM50(opt.model_name, job.init_soc, {'CC', job.crate}, ...
                common_args{:}, ...
                'VEnd', job.v_end);

        case 'sequence'
            % Example for future extension:
            % job.profile can be {'REST',0.05,'FUDS','REST',0.05,'BJDST'}
            sim = sim_profile_sequence_LGM50(opt.model_name, job.init_soc, job.profile, ...
                common_args{:}, ...
                'VMin', opt.VMin, 'VMax', opt.VMax, ...
                'ShowProgress', false);

        otherwise
            error('Unknown job.kind: %s', job.kind);
    end
end

function data = convert_sim_to_pinn_dataset(sim, opt)
    p = sim.p;
    NT = numel(sim.time);

    data = struct();
    data.time = sim.time(:);
    data.current = sim.current(:);
    data.voltage = sim.voltage(:);
    data.soc = sim.soc(:);

    % States from sim_profile_LGM50 are stored as [state x time].
    data.c_s_n = sim.c_s_n;       % (Nr-1) x NT
    data.c_s_p = sim.c_s_p;       % (Nr-1) x NT
    data.c_e = sim.c_e;           % (Nx-3) x NT

    data.css_n = data.c_s_n(end, :).';
    data.css_p = data.c_s_p(end, :).';

    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cs_ave_n')
        data.cs_ave_n = sim.spme_outputs.cs_ave_n(:);
        data.cs_ave_p = sim.spme_outputs.cs_ave_p(:);
    else
        data.cs_ave_n = radial_average_spherical(data.c_s_n, p);
        data.cs_ave_p = radial_average_spherical(data.c_s_p, p);
    end

    % Electrolyte boundary/interface values and region means, if available.
    data.ce0n = nan(NT,1); data.cens = nan(NT,1);
    data.cesp = nan(NT,1); data.ce0p = nan(NT,1);
    data.cen_bar = nan(NT,1); data.ces_bar = nan(NT,1); data.cep_bar = nan(NT,1);
    if isfield(sim, 'spme_outputs')
        f = sim.spme_outputs;
        names = {'ce0n','cens','cesp','ce0p','cen_bar','ces_bar','cep_bar'};
        for k = 1:numel(names)
            if isfield(f, names{k})
                data.(names{k}) = f.(names{k})(:);
            end
        end
    end

    % State order conversion for PINN-SPMe.
    % MATLAB: [c_s_n; c_s_p; c_e]
    % PINN:   [c_s_n; c_e; c_s_p]
    data.x_matlab_order = sim.x;
    data.x_pinn_order = [data.c_s_n; data.c_e; data.c_s_p];

    % Core matrix similar to old SPM_FDM CSV style:
    % [t, I, V, css_n, cs_ave_n, css_p, cs_ave_p, c_s_n, c_e, c_s_p]
    data.matrix = [ ...
        data.time, ...
        data.current, ...
        data.voltage, ...
        data.soc, ...
        data.css_n, ...
        data.cs_ave_n, ...
        data.css_p, ...
        data.cs_ave_p, ...
        data.ce0n, data.cens, data.cesp, data.ce0p, ...
        data.cen_bar, data.ces_bar, data.cep_bar, ...
        data.c_s_n.', ...
        data.c_e.', ...
        data.c_s_p.' ...
    ];

    data.column_names = build_column_names(opt, p);
    data.discretization = struct('Nr', p.Nr, 'Nxn', p.Nxn, 'Nxs', p.Nxs, ...
                                 'Nxp', p.Nxp, 'Nx', p.Nx);
    data.sign_convention = '+ current = discharge, - current = charge';
end

function avg = radial_average_spherical(c_s, p)
    % c_s is (Nr-1) x NT for interior radial nodes only.
    % This is only fallback; prefer SPMe optional output cs_ave_n/p.
    Nr = p.Nr;
    dr = 1 / Nr;
    r = (dr:dr:(1-dr)).';
    w = r.^2;
    avg = (w.' * c_s ./ sum(w)).';
end

function names = build_column_names(opt, p)
    names = ["time_s","current_A","voltage_V","soc", ...
             "css_n","cs_ave_n","css_p","cs_ave_p", ...
             "ce0n","cens","cesp","ce0p", ...
             "cen_bar","ces_bar","cep_bar"];
    for i = 1:(p.Nr-1)
        names(end+1) = "c_s_n_" + i; %#ok<AGROW>
    end
    for i = 1:(p.Nx-3)
        names(end+1) = "c_e_" + i; %#ok<AGROW>
    end
    for i = 1:(p.Nr-1)
        names(end+1) = "c_s_p_" + i; %#ok<AGROW>
    end
end

function write_core_csv(path, data)
    T = array2table(data.matrix, 'VariableNames', cellstr(data.column_names));
    writetable(T, path);
end

function write_summary_csv(path, summary)
    if isempty(summary)
        return;
    end
    T = struct2table(summary);
    writetable(T, path);
end

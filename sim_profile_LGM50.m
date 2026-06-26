function sim = sim_profile_LGM50(model_name, init_soc, profile_def, varargin)
%SIM_PROFILE_LGM50 Simulate DFN, SPMe, or ECM on CALCE-derived profiles.
%   sim = sim_profile_LGM50(model_name, init_soc, profile_def, ...)
%
%   Inputs:
%     model_name   : 'DFN', 'SPMe', or 'ECM' (case-insensitive)
%     init_soc     : initial SoC in [0,1], default 0.9
%     profile_def  : single profile definition, e.g.
%                    'DST'
%                    {'CC', -1}
%                    {'REST', 0.1}  % rest for 0.1 hour at I = 0 A
%
%   Name-Value options:
%     'Nxn'        : # x-nodes in anode electrode (default 70)
%     'Nxs'        : # x-nodes in separator (default 35)
%     'Nxp'        : # x-nodes in cathode electrode (default 70)
%     'Nr'         : # particle-state nodes for SPMe/DFN init (default 30)
%     'DeltaT'     : simulation timestep [s] (default 1)
%     'Tamb'       : ambient temperature [degC] (default 25)
%     'SaveOutput' : true/false, save .mat to results folder (default false)
%     'OutputFile' : optional full/relative .mat path when SaveOutput=true
%     'ECMOrder'   : number of RC branches for ECM (default 2)
%     'Crate'      : C-rate for CC mode (default 1).
%                    sign convention: + discharge, - charge
%     'DegradationModel' : SPMe degradation model
%                    {'none','ramadass','yang','safari'} (default 'none')
%     'ElectrolyteModel' : SPMe electrolyte dynamics model
%                    {'fdm','fvm'} (default 'fdm')
%     'UseElectrolyteActivity' : include electrolyte activity correction
%                    factor in SPMe electrolyte-polarization voltage term
%                    (default false)
%     'EnableLithiumPlating' : enable lithium-plating side reaction add-on
%                    for SEI models (default false)
%     'EnableAnodePorosityLoss' : enable anode porosity-loss add-on from
%                    film growth for SEI models (default false)
%     'DFNLplCouplingMode' : DFN Li-plating film-drop coupling mode,
%                    {'spme_compatible','yang_total_flux'}
%                    (default 'spme_compatible')
%     'Nsei'       : # interior EC concentration nodes in SEI for Safari
%                    model only (default 10)
%     'EpsilonSEI' : SEI/electrolyte EC partition factor for Safari
%                    model only (default 1.0)
%     'InitialState' : state struct for continuation from prior segment.
%                    For SPMe: struct('x', x_last)
%                    For DFN : struct('x', x_last, 'z', z_last)
%                    For ECM : struct('x', x_last)
%     'SuppressSummary' : true/false to suppress console summary print
%                    (default false)
%     'DiagnosticsFile' : optional .mat path for failure diagnostics
%     'DiagnosticContext' : optional struct copied into diagnostics
%     'SPMeParamOverrides' : optional struct of SPMe parameter overrides
%     'VEnd'       : custom terminal voltage [V]. If provided:
%                    - dynamic profiles and CC discharge stop when V <= VEnd
%                    - CC charge stops when V >= VEnd
%                    - dynamic profiles repeat until VEnd is reached
%
%   Example: custom discretization + save
%   sim = sim_profile_LGM50('SPMe',0.9,'FUDS', ...
%       'Nxn', 80, 'Nxs', 40, 'Nxp', 80, 'Nr', 35, ...
%       'SaveOutput', true);
%
%   Notes:
%   - Uses profile files in data_CALCE: <PROFILE>_LGM50.csv
%   - Current sign is used as-is from the CSV.

    if nargin < 1 || isempty(model_name)
        model_name = 'SPMe';
    end
    if nargin < 2 || isempty(init_soc)
        init_soc = 0.9;
    end
    if nargin < 3 || isempty(profile_def)
        profile_def = 'DST';
    end

    model_name = upper(char(model_name));
    if ~ismember(model_name, {'DFN', 'SPME', 'ECM'})
        error('model_name must be ''DFN'', ''SPMe'', or ''ECM''.');
    end

    parser = inputParser;
    addParameter(parser, 'Nxn', 20);
    addParameter(parser, 'Nxs', 10);
    addParameter(parser, 'Nxp', 20);
    addParameter(parser, 'Nr', 20);
    addParameter(parser, 'DeltaT', 1);
    addParameter(parser, 'Tamb', 25);
    addParameter(parser, 'SaveOutput', false);
    addParameter(parser, 'OutputFile', '');
    addParameter(parser, 'ECMOrder', 2);
    addParameter(parser, 'Crate', 1);
    addParameter(parser, 'DegradationModel', 'none');
    addParameter(parser, 'ElectrolyteModel', 'fdm');
    addParameter(parser, 'UseElectrolyteActivity', false);
    addParameter(parser, 'EnableLithiumPlating', []);
    addParameter(parser, 'EnableAnodePorosityLoss', []);
    addParameter(parser, 'DFNLplCouplingMode', 'spme_compatible');
    addParameter(parser, 'Nsei', 10); % added for Safari's SEI model
    addParameter(parser, 'EpsilonSEI', 0.05);  % added for Safari's SEI model 
    addParameter(parser, 'VEnd', []);
    addParameter(parser, 'InitialState', struct());
    addParameter(parser, 'SuppressSummary', false);
    addParameter(parser, 'SPMeOutputs', {});
    addParameter(parser, 'DiagnosticsFile', '');
    addParameter(parser, 'DiagnosticContext', struct());
    addParameter(parser, 'SPMeParamOverrides', struct());
    parse(parser, varargin{:});
    opt = parser.Results;
    has_initial_state = isstruct(opt.InitialState) && ...
        isfield(opt.InitialState, 'x') && ~isempty(opt.InitialState.x);
    has_initial_params = has_initial_state && ...
        isfield(opt.InitialState, 'p') && ~isempty(opt.InitialState.p);
    use_spme_continuation_params = strcmp(model_name, 'SPME') && has_initial_params;

    if ~use_spme_continuation_params && ...
            (~isscalar(init_soc) || init_soc < 0 || init_soc > 1)
        error('init_soc must be a scalar in [0, 1].');
    end

    base_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(base_dir, 'param'));
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'DFN'));
    ecm_dir = fullfile(base_dir, 'ECM');
    if isfolder(ecm_dir)
        addpath(ecm_dir);
    end
    ecm2_dir = fullfile(base_dir, 'ECM 2');
    if isfolder(ecm2_dir)
        addpath(ecm2_dir);
    end

    if ~isscalar(opt.Crate) || ~isfinite(opt.Crate)
        error('Crate must be a finite scalar.');
    end
    degradation_model = lower(char(opt.DegradationModel));
    if ~ismember(degradation_model, {'none', 'ramadass', 'yang', 'safari'})
        error('DegradationModel must be one of: none, ramadass, yang, safari.');
    end
    electrolyte_model = lower(char(opt.ElectrolyteModel));
    if ~ismember(electrolyte_model, {'fdm', 'fvm'})
        error('ElectrolyteModel must be one of: fdm, fvm.');
    end
    if ~isscalar(opt.Nsei) || ~isfinite(opt.Nsei) || opt.Nsei < 1 || abs(opt.Nsei - round(opt.Nsei)) > 0  % added for Safari's SEI model
        error('Nsei must be a positive integer.');
    end
    if ~isscalar(opt.EpsilonSEI) || ~isfinite(opt.EpsilonSEI) || opt.EpsilonSEI < 0   % added for Safari's SEI model
        error('EpsilonSEI must be a nonnegative finite scalar.');
    end
    if ~(isscalar(opt.UseElectrolyteActivity) || islogical(opt.UseElectrolyteActivity))
        error('UseElectrolyteActivity must be a scalar logical/numeric value.');
    end
    use_electrolyte_activity = logical(opt.UseElectrolyteActivity);
    if isempty(opt.EnableLithiumPlating)
        enable_lithium_plating = false;
    elseif isscalar(opt.EnableLithiumPlating) && (isnumeric(opt.EnableLithiumPlating) || islogical(opt.EnableLithiumPlating))
        enable_lithium_plating = logical(opt.EnableLithiumPlating);
    else
        error('EnableLithiumPlating must be empty or a scalar logical/numeric value.');
    end
    if isempty(opt.EnableAnodePorosityLoss)
        enable_anode_porosity_loss = false;
    elseif isscalar(opt.EnableAnodePorosityLoss) && (isnumeric(opt.EnableAnodePorosityLoss) || islogical(opt.EnableAnodePorosityLoss))
        enable_anode_porosity_loss = logical(opt.EnableAnodePorosityLoss);
    else
        error('EnableAnodePorosityLoss must be empty or a scalar logical/numeric value.');
    end
    if strcmp(degradation_model, 'none')
        enable_lithium_plating = false;
        enable_anode_porosity_loss = false;
    end
    dfn_lpl_coupling_mode = lower(char(opt.DFNLplCouplingMode));
    if ~ismember(dfn_lpl_coupling_mode, {'spme_compatible', 'yang_total_flux'})
        error('DFNLplCouplingMode must be one of: spme_compatible, yang_total_flux.');
    end
    use_v_end = ~isempty(opt.VEnd);
    if use_v_end && (~isscalar(opt.VEnd) || ~isfinite(opt.VEnd))
        error('VEnd must be a finite scalar when provided.');
    end

    [profile_name, crate_single, rest_hours_single] = parse_single_profile(profile_def, opt.Crate);
    if strcmp(profile_name, 'CC')
        opt.Crate = crate_single;
    end

    is_cc_profile = strcmp(profile_name, 'CC');
    is_rest_profile = strcmp(profile_name, 'REST');
    is_cc_charge = is_cc_profile && (opt.Crate < 0);
    is_dynamic_profile = ~(is_cc_profile || is_rest_profile);
    repeat_dynamic_for_v_end = is_dynamic_profile && use_v_end;
    max_dynamic_repeats = 200;
    if is_cc_profile
        if opt.Crate == 0
            error('For CC mode, Crate must be non-zero.');
        end
        % For CC modes, generate a long-enough horizon and terminate by voltage.
        max_time_s = ceil(1.5 * 3600 / abs(opt.Crate));
        NT = max(2, floor(max_time_s / max(opt.DeltaT, eps)) + 1);
        t_profile = (0:NT-1)';
        t_profile_full = t_profile;
        I = [];
    elseif is_rest_profile
        rest_duration_s = rest_hours_single * 3600;
        NT = max(2, floor(rest_duration_s / max(opt.DeltaT, eps)) + 1);
        t_profile = (0:NT-1)' * opt.DeltaT;
        t_profile_full = t_profile;
        I = zeros(NT, 1);
    else
        profile_path = fullfile(base_dir, 'data_CALCE', [profile_name '_LGM50.csv']);
        if ~isfile(profile_path)
            error('Profile file not found: %s', profile_path);
        end

        opts = detectImportOptions(profile_path, 'VariableNamingRule', 'preserve');
        prof = readtable(profile_path, opts);
        t_profile_base = prof.('Test_Time(s)');
        I_base = prof.('Current(A)');
        I_base = I_base(:);
        NT_base = numel(I_base);
        if NT_base < 2
            error('Profile file must contain at least 2 rows.');
        end
        if repeat_dynamic_for_v_end
            NT = NT_base * max_dynamic_repeats;
            I = repmat(I_base, max_dynamic_repeats, 1);
            t_profile_full = repmat(t_profile_base, max_dynamic_repeats, 1);
        else
            NT = NT_base;
            I = I_base;
            t_profile_full = t_profile_base;
        end
        t_profile = t_profile_full;
    end

    % Fixed-step simulation settings.
    delta_t = opt.DeltaT; % [s]
    if ~isscalar(delta_t) || delta_t <= 0
        error('DeltaT must be a positive scalar.');
    end
    t = (0:NT-1)' * delta_t;
    t_amb = opt.Tamb; % [degC]
    disable_thermal = true;

    Nr = opt.Nr;
    Nxn = opt.Nxn;
    Nxs = opt.Nxs;
    Nxp = opt.Nxp;
    if any([Nr, Nxn, Nxs, Nxp] <= 2)
        error('Nr, Nxn, Nxs, and Nxp must be > 2.');
    end

    % Requested initialization sequence:
    % SoC -> Veq by interpolation; then Veq -> (csn0, csp0) via init_cs.
    % Continuation with InitialState.p skips fresh-cell SoC initialization.
    if use_spme_continuation_params
        calc = struct('Veq', NaN);
    else
        calc = calculate_equilibrium_voltage_from_soc(init_soc);
    end

    switch model_name
        case 'SPME'
            if use_spme_continuation_params
                p = opt.InitialState.p;
                p.delta_t = delta_t;
                p.T_amb = t_amb + 273.15;
                p.disable_thermal = disable_thermal;
                p.degradation_model = degradation_model;
                p.electrolyte_model = electrolyte_model;
                p.use_electrolyte_activity = use_electrolyte_activity;
                p.enable_lithium_plating = enable_lithium_plating;
                p.enable_anode_porosity_loss = enable_anode_porosity_loss;
                Nr = p.Nr;
                Nxn = p.Nxn;
                Nxs = p.Nxs;
                Nxp = p.Nxp;
            else
                p = load_SPMe_params(init_soc, t_amb, Nr, delta_t, Nxn, Nxs, Nxp, disable_thermal, degradation_model, electrolyte_model, use_electrolyte_activity, enable_lithium_plating, enable_anode_porosity_loss);
                if strcmpi(degradation_model, 'safari') % added for Safari's SEI model
                    p.Nsei = round(opt.Nsei);
                    p.epsilon_sei = opt.EpsilonSEI;
                end
            end
            p = apply_spme_param_overrides(p, opt.SPMeParamOverrides);
            if is_cc_profile
                I = (p.OneC * opt.Crate) * ones(NT, 1);
            end

            if use_spme_continuation_params
                x0 = opt.InitialState.x(:);
                csn0 = NaN;
                csp0 = NaN;
            else
                [csn0, csp0] = init_cs(p, calc.Veq);
                p.csn0 = csn0;
                p.csp0 = csp0;

                c_n0 = csn0 * ones(p.Nr - 1, 1);
                c_p0 = csp0 * ones(p.Nr - 1, 1);
                c_e0 = p.c_e * ones(p.Nx - 3, 1);
                x0 = [c_n0; c_p0; c_e0];
                if ~strcmpi(char(p.degradation_model), 'none')
                    delta_sei0 = 0;
                    x0 = [x0; delta_sei0];
                end
                if p.enable_lithium_plating
                    delta_li0 = 0;
                    x0 = [x0; delta_li0];
                end
                if p.enable_anode_porosity_loss
                    epsilon_e_n0 = p.epsilon_e_n;
                    x0 = [x0; epsilon_e_n0];
                end
                if strcmpi(char(p.degradation_model), 'safari') % addd for Safari's model
                    c_ec_sei0 = p.epsilon_sei * p.c0_ec * ones(p.Nsei, 1);
                    x0 = [x0; c_ec_sei0];
                end
                nx_expected = numel(x0);
                if has_initial_state
                    x0 = opt.InitialState.x(:);
                    if numel(x0) ~= nx_expected
                        error('InitialState.x size mismatch for SPMe continuation. Expected %d states, got %d.', nx_expected, numel(x0));
                    end
                end
            end
            p.x0 = x0;

            spme_solver = 'idas';
            if strcmpi(char(p.degradation_model), 'safari')
                % Safari's added SEI diffusion states are more robust with cvodes.
                spme_solver = 'cvodes';
            end
            spme_ss = SPMe(p, spme_solver);
            if use_spme_continuation_params && numel(x0) ~= numel(spme_ss.x)
                error('InitialState.x size mismatch for SPMe continuation. Expected %d states, got %d.', numel(spme_ss.x), numel(x0));
            end
            req_spme_out = normalize_name_list(opt.SPMeOutputs);
            spme_out_vals = [];
            if ~isempty(req_spme_out)
                spme_out_vals = zeros(numel(req_spme_out), NT);
                out0 = spme_ss.optional_outputs(x0, I(1), req_spme_out);
                for ii = 1:numel(req_spme_out)
                    spme_out_vals(ii, 1) = out0.(req_spme_out{ii});
                end
            end

            x = zeros(length(x0), NT);
            y0 = spme_ss.h('x', x0, 'p', I(1));
            y = zeros(length(full(y0.y)), NT);
            x(:, 1) = x0;
            y(:, 1) = full(y0.y);
            t_last = NT;
            terminated_early = false;
            termination_reason = '';

            for k = 1:NT-1
                try
                    [x(:, k+1), y(:, k+1)] = spme_ss.step(x(:, k), I(k));
                catch ME
                    diagnostics_file = save_spme_step_failure_diagnostics( ...
                        opt.DiagnosticsFile, opt.DiagnosticContext, ME, ...
                        profile_name, init_soc, p, spme_solver, t, t_profile, I, ...
                        x, y, k, req_spme_out, spme_out_vals);
                    if ~isempty(diagnostics_file)
                        fprintf(2, ['SPMe step diagnostics saved before rethrowing solver failure:\n' ...
                            '  %s\n'], diagnostics_file);
                    end
                    rethrow(ME);
                end
                if ~isempty(req_spme_out)
                    outk = spme_ss.optional_outputs(x(:, k+1), I(k), req_spme_out);
                    for ii = 1:numel(req_spme_out)
                        spme_out_vals(ii, k+1) = outk.(req_spme_out{ii});
                    end
                end
                Vk = y(1, k+1);
                if use_v_end
                    if is_cc_charge
                        if Vk >= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for charge (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    else
                        if Vk <= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for discharge/dynamic (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    end
                else
                    if Vk <= p.volt_min
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmin (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_min, t(t_last), t_last);
                        break;
                    elseif Vk >= p.volt_max
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmax (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_max, t(t_last), t_last);
                        break;
                    end
                end
            end
            if use_v_end && ~terminated_early
                termination_reason = sprintf('Profile horizon ended before reaching VEnd=%.6f.', opt.VEnd);
            end

            sim = struct();
            sim.model = model_name;
            sim.profile = profile_name;
            sim.profile_time = t_profile(1:t_last);
            sim.time = t(1:t_last);
            sim.current = I(1:t_last);
            sim.voltage = y(1, 1:t_last).';
            sim.soc = y(2, 1:t_last).';
            sim.c_s_n = x(1:p.Nr-1, 1:t_last);
            sim.c_s_p = x(p.Nr:2*(p.Nr-1), 1:t_last);
            idx_ce_start = 2*p.Nr - 1;
            idx_ce_end = idx_ce_start + p.Nx - 4;
            sim.c_e = x(idx_ce_start:idx_ce_end, 1:t_last);
            idx_next = idx_ce_end + 1;
            if ~strcmpi(char(p.degradation_model), 'none')
                idx_delta = idx_next;
                idx_next = idx_next + 1;
                sim.delta_sei = x(idx_delta, 1:t_last).';
            end
            if p.enable_lithium_plating
                idx_delta_li = idx_next;
                idx_next = idx_next + 1;
                sim.delta_li = x(idx_delta_li, 1:t_last).';
            end
            if p.enable_anode_porosity_loss
                idx_epsilon_e_n = idx_next;
                idx_next = idx_next + 1;
                sim.epsilon_e_n = x(idx_epsilon_e_n, 1:t_last).';
            end
            if strcmpi(char(p.degradation_model), 'safari')
                idx_c_ec_start = idx_next;
                idx_c_ec_end = idx_c_ec_start + p.Nsei - 1;
                sim.c_ec_sei = x(idx_c_ec_start:idx_c_ec_end, 1:t_last);
                safari_state_indices = struct('c_ec_sei', idx_c_ec_start:idx_c_ec_end);
                if exist('idx_delta', 'var')
                    safari_state_indices.delta_sei = idx_delta;
                end
                sim.safari = struct('Nsei', p.Nsei, 'epsilon_sei', p.epsilon_sei, ...
                                    'state_indices', safari_state_indices);
            end
            sim.x = x(:, 1:t_last);
            sim.y = y(:, 1:t_last);
            if ~isempty(req_spme_out)
                sim.spme_outputs = struct();
                for ii = 1:numel(req_spme_out)
                    sim.spme_outputs.(req_spme_out{ii}) = spme_out_vals(ii, 1:t_last).';
                end
            end
            sim.init_soc_input = init_soc;
            sim.init_voltage_from_soc = calc.Veq;
            sim.init_csn0 = csn0;
            sim.init_csp0 = csp0;
            sim.p = p;
            sim.discretization = struct('Nr', p.Nr, 'Nxn', p.Nxn, 'Nxs', p.Nxs, 'Nxp', p.Nxp);
            sim.terminated_early = terminated_early;
            sim.termination_reason = termination_reason;
            sim.t_last = t_last;

        case 'DFN'
            if ~ismember(degradation_model, {'none', 'safari'})
                error('DFN degradation currently supports DegradationModel values: none, safari.');
            end
            p = load_DFN_params(init_soc, t_amb, delta_t, disable_thermal, ...
                Nr, Nxn, Nxs, Nxp, degradation_model, enable_lithium_plating, ...
                enable_anode_porosity_loss, dfn_lpl_coupling_mode);
            if is_cc_profile
                I = (p.OneC * opt.Crate) * ones(NT, 1);
            end

            [csn0, csp0] = init_cs(p, calc.Veq);
            p.csn0 = csn0;
            p.csp0 = csp0;

            % Rebuild DFN initial states using concentrations from init_cs.
            c_s_n0 = zeros(p.PadeOrder, 1);
            c_s_p0 = zeros(p.PadeOrder, 1);
            c_s_n0(3) = csn0;
            c_s_p0(3) = csp0;

            c_s_n = repmat(c_s_n0, [p.Nxn-1, 1]);
            c_s_p = repmat(c_s_p0, [p.Nxp-1, 1]);
            c_e = p.c_e * ones(p.Nx-3, 1);
            T = t_amb + 273.15;

            Uref_n0 = refPotentialAnode(p, csn0 / p.c_s_n_max) * ones(p.Nxn-1, 1);
            Uref_p0 = refPotentialCathode(p, csp0 / p.c_s_p_max) * ones(p.Nxp-1, 1);
            i_en = zeros(p.Nxn-1, 1);
            i_ep = zeros(p.Nxp-1, 1);
            phi_e = zeros(p.Nx-1, 1);
            jn = zeros(p.Nxn-1, 1);
            jp = zeros(p.Nxp-1, 1);

            if isfield(p, 'dfn_degradation_enabled') && p.dfn_degradation_enabled
                delta_sei_n0 = zeros(p.Nxn-1, 1);
                delta_li_n0 = zeros(p.Nxn-1, 1);
                epsilon_e_n_dyn0 = p.epsilon_e_n * ones(p.Nxn-1, 1);
                p.x0 = [c_s_n; c_s_p; c_e; delta_sei_n0; delta_li_n0; epsilon_e_n_dyn0; T];
            else
                p.x0 = [c_s_n; c_s_p; c_e; T];
            end
            p.z0 = [Uref_n0; Uref_p0; i_en; i_ep; phi_e; jn; jp];
            [~, ~, p.y0] = dae_dfn(p.x0, p.z0, 0, p);
            expected_dfn_x_len = numel(p.x0);
            if isstruct(opt.InitialState) && ~isempty(opt.InitialState)
                if isfield(opt.InitialState, 'x') && ~isempty(opt.InitialState.x)
                    p.x0 = opt.InitialState.x(:);
                    if numel(p.x0) ~= expected_dfn_x_len
                        error('InitialState.x size mismatch for DFN continuation.');
                    end
                end
                if isfield(opt.InitialState, 'z') && ~isempty(opt.InitialState.z)
                    p.z0 = opt.InitialState.z(:);
                    if numel(p.z0) ~= numel([Uref_n0; Uref_p0; i_en; i_ep; phi_e; jn; jp])
                        error('InitialState.z size mismatch for DFN continuation.');
                    end
                end
                [~, ~, p.y0] = dae_dfn(p.x0, p.z0, 0, p);
            end

            x = zeros(length(p.x0), NT);
            z = zeros(length(p.z0), NT);
            y = zeros(length(p.y0), NT);
            x(:, 1) = p.x0;
            z(:, 1) = p.z0;
            y(:, 1) = p.y0;
            t_last = NT;
            terminated_early = false;
            termination_reason = '';

            for k = 1:NT-1
                if k == 1
                    Cur_vec = [I(k), I(k), I(k+1)];
                else
                    Cur_vec = [I(k-1), I(k), I(k+1)];
                end
                [x(:, k+1), z(:, k+1), y(:, k+1)] = step_dfn(x(:, k), z(:, k), Cur_vec, p, 0);
                Vk = y(1, k+1);
                if use_v_end
                    if is_cc_charge
                        if Vk >= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for charge (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    else
                        if Vk <= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for discharge/dynamic (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    end
                else
                    if Vk <= p.volt_min
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmin (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_min, t(t_last), t_last);
                        break;
                    elseif Vk >= p.volt_max
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmax (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_max, t(t_last), t_last);
                        break;
                    end
                end
            end
            if use_v_end && ~terminated_early
                termination_reason = sprintf('Profile horizon ended before reaching VEnd=%.6f.', opt.VEnd);
            end

            [c_s_n_flat, c_s_p_flat, c_e_hist, T_hist, dfn_xdeg] = parse_x(x, t_last, p);
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, dfn_ydeg] = parse_y(y, t_last, p);
            c_s_n_pade = zeros(p.PadeOrder, p.Nxn-1, t_last);
            c_s_p_pade = zeros(p.PadeOrder, p.Nxp-1, t_last);
            for k = 1:t_last
                c_s_n_pade(:, :, k) = reshape(c_s_n_flat(:, k), p.PadeOrder, p.Nxn-1);
                c_s_p_pade(:, :, k) = reshape(c_s_p_flat(:, k), p.PadeOrder, p.Nxp-1);
            end

            sim = struct();
            sim.model = model_name;
            sim.profile = profile_name;
            sim.profile_time = t_profile(1:t_last);
            sim.time = t(1:t_last);
            sim.current = I(1:t_last);
            sim.voltage = y(1, 1:t_last).';
            sim.soc = y(2, 1:t_last).';
            sim.c_s_n_pade = c_s_n_pade;
            sim.c_s_p_pade = c_s_p_pade;
            sim.c_e = c_e_hist;
            sim.temperature = T_hist;
            sim.x = x(:, 1:t_last);
            sim.z = z(:, 1:t_last);
            sim.y = y(:, 1:t_last);
            sim.dfn_outputs = dfn_ydeg;
            sim.delta_sei_n = dfn_xdeg.delta_sei_n.';
            sim.delta_li_n = dfn_xdeg.delta_li_n.';
            sim.epsilon_e_n_dyn = dfn_xdeg.epsilon_e_n_dyn.';
            sim.dfn_lpl_coupling_mode = p.dfn_lpl_coupling_mode;
            sim.init_soc_input = init_soc;
            sim.init_voltage_from_soc = calc.Veq;
            sim.init_csn0 = csn0;
            sim.init_csp0 = csp0;
            sim.p = p;
            sim.discretization = struct('PadeOrder', p.PadeOrder, 'Nr', p.Nr, 'Nxn', p.Nxn, 'Nxs', p.Nxs, 'Nxp', p.Nxp);
            sim.terminated_early = terminated_early;
            sim.termination_reason = termination_reason;
            sim.t_last = t_last;

        case 'ECM'
            n_order = opt.ECMOrder;
            if ~isscalar(n_order) || n_order < 1 || abs(n_order - round(n_order)) > 0
                error('ECMOrder must be a positive integer.');
            end

            [p, ecm_p, ocv_p] = load_ECM_params(init_soc, n_order);
            p.delta_t = delta_t;
            p.T_amb = t_amb + 273.15;
            if is_cc_profile
                I = (p.OneC * opt.Crate) * ones(NT, 1);
            end

            [csn0, csp0] = init_cs(p, calc.Veq);
            [cn_low, cp_low] = init_cs(p, p.volt_min);
            [cn_high, cp_high] = init_cs(p, p.volt_max);

            z0_n = (csn0 - cn_low) / max(cn_high - cn_low, eps);
            z0_p = (csp0 - cp_low) / max(cp_high - cp_low, eps);
            z0_from_initcs = min(max(0.5 * (z0_n + z0_p), 0), 1);

            % Match legacy ECM behavior: initialize SOC state directly from
            % requested init_soc (as in sim_ECM_CC / load_ECM_params).
            z0 = min(max(init_soc, 0), 1);
            x0 = p.x0;
            x0(1) = z0;
            if isstruct(opt.InitialState) && isfield(opt.InitialState, 'x') && ~isempty(opt.InitialState.x)
                x0 = opt.InitialState.x(:);
                if numel(x0) ~= numel(p.x0)
                    error('InitialState.x size mismatch for ECM continuation.');
                end
            end
            ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');

            x = zeros(length(x0), NT);
            V = zeros(1, NT);
            x(:, 1) = x0;
            V(1, 1) = p.y0;
            t_last = NT;
            terminated_early = false;
            termination_reason = '';

            for k = 1:NT-1
                [x(:, k+1), V(:, k+1)] = ecm_ss.step(x(:, k), I(k));
                Vk = V(1, k+1);
                if use_v_end
                    if is_cc_charge
                        if Vk >= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for charge (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    else
                        if Vk <= opt.VEnd
                            t_last = k+1;
                            terminated_early = true;
                            termination_reason = sprintf('Early termination: voltage reached VEnd for discharge/dynamic (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                                Vk, opt.VEnd, t(t_last), t_last);
                            break;
                        end
                    end
                else
                    if Vk <= p.volt_min
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmin (%.6f <= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_min, t(t_last), t_last);
                        break;
                    elseif Vk >= p.volt_max
                        t_last = k+1;
                        terminated_early = true;
                        termination_reason = sprintf('Early termination: voltage hit Vmax (%.6f >= %.6f) at t=%.2f s, step=%d.', ...
                            Vk, p.volt_max, t(t_last), t_last);
                        break;
                    end
                end
            end
            if use_v_end && ~terminated_early
                termination_reason = sprintf('Profile horizon ended before reaching VEnd=%.6f.', opt.VEnd);
            end

            sim = struct();
            sim.model = model_name;
            sim.profile = profile_name;
            sim.profile_time = t_profile(1:t_last);
            sim.time = t(1:t_last);
            sim.current = I(1:t_last);
            sim.voltage = V(1, 1:t_last).';
            sim.soc = x(1, 1:t_last).';
            sim.x = x(:, 1:t_last);
            sim.init_soc_input = init_soc;
            sim.init_voltage_from_soc = calc.Veq;
            sim.init_csn0 = csn0;
            sim.init_csp0 = csp0;
            sim.init_soc_used = z0;
            sim.init_soc_from_init_cs = z0_from_initcs;
            sim.init_soc_from_init_cs_anode = z0_n;
            sim.init_soc_from_init_cs_cathode = z0_p;
            sim.ocv_p = ocv_p;
            sim.ecm_p = ecm_p;
            sim.p = p;
            sim.discretization = struct('ECMOrder', n_order);
            sim.terminated_early = terminated_early;
            sim.termination_reason = termination_reason;
            sim.t_last = t_last;
    end

    if opt.SaveOutput
        if isempty(opt.OutputFile)
            out_dir = fullfile(base_dir, 'results', model_name);
            if ~isfolder(out_dir)
                mkdir(out_dir);
            end
            out_name = sprintf('%s_%s_%0.3fSOC.mat', lower(model_name), profile_name, init_soc);
            out_path = fullfile(out_dir, out_name);
        else
            out_path = opt.OutputFile;
        end
        save(out_path, 'sim');
        sim.saved_file = out_path;
    else
        sim.saved_file = '';
    end

    if ~opt.SuppressSummary
        fprintf('\nSim complete: model=%s, profile=%s, init SoC=%.3f\n', model_name, profile_name, init_soc);
        fprintf('Initial voltage from SoC interpolation: %.6f V\n', calc.Veq);
        fprintf('Initial concentrations from init_cs(Veq): csn0=%.6f, csp0=%.6f [mol/m^3]\n', sim.init_csn0, sim.init_csp0);
        if strcmp(model_name, 'ECM')
            fprintf('ECM order (RC branches): %d\n', opt.ECMOrder);
        else
            fprintf('Discretization (Nxn,Nxs,Nxp,Nr): (%d,%d,%d,%d)\n', Nxn, Nxs, Nxp, Nr);
            if strcmpi(degradation_model, 'safari')
                fprintf('Safari SEI discretization: Nsei=%d, EpsilonSEI=%.6g\n', round(opt.Nsei), opt.EpsilonSEI);
            end
            if enable_lithium_plating || enable_anode_porosity_loss
                fprintf('SEI add-ons: LithiumPlating=%d, AnodePorosityLoss=%d\n', ...
                    enable_lithium_plating, enable_anode_porosity_loss);
            end
        end
        if sim.terminated_early
            fprintf('%s\n', sim.termination_reason);
        else
            fprintf('Completed full profile without voltage-limit termination.\n');
        end
        if opt.SaveOutput
            fprintf('Saved output file: %s\n', sim.saved_file);
        else
            fprintf('Output not saved to disk (SaveOutput=false).\n');
        end
        fprintf('Final voltage: %.6f V\n\n', sim.voltage(end));
    end
end

function [profile_name, crate_single, rest_hours_single] = parse_single_profile(profile_def, default_crate)
    valid_dynamic = {'BJDST', 'DST', 'FUDS', 'US06'};
    crate_single = NaN;
    rest_hours_single = NaN;
    if ischar(profile_def) || isstring(profile_def)
        profile_name = upper(char(profile_def));
        if strcmp(profile_name, 'CC')
            if ~isscalar(default_crate) || ~isfinite(default_crate) || default_crate == 0
                error('For profile_def=''CC'', Crate must be provided as a non-zero scalar.');
            end
            crate_single = default_crate;
            return;
        end
        if strcmp(profile_name, 'REST')
            error('For REST mode, use profile_def={''REST'', duration_hours}.');
        end
        if ~ismember(profile_name, valid_dynamic)
            error('profile_def must be one of BJDST, DST, FUDS, US06, ''CC'', {''CC'', crate}, or {''REST'', duration_hours}.');
        end
        return;
    end

    if ~iscell(profile_def)
        error('profile_def must be char/string for one profile, or {''CC'', crate}/{''REST'', hours}.');
    end

    if numel(profile_def) ~= 2 || ~(ischar(profile_def{1}) || isstring(profile_def{1})) || ...
            ~isnumeric(profile_def{2}) || ~isscalar(profile_def{2}) || ...
            ~isfinite(profile_def{2})
        error(['sim_profile_LGM50 supports a single profile only. ', ...
               'Use sim_profile_sequence_LGM50 for concatenated profiles, ', ...
               'or pass profile_def={''CC'',crate} / {''REST'',hours}.']);
    end
    token = upper(char(profile_def{1}));
    if strcmp(token, 'CC')
        if profile_def{2} == 0
            error('For {''CC'', crate}, crate must be non-zero.');
        end
        profile_name = 'CC';
        crate_single = profile_def{2};
    elseif strcmp(token, 'REST')
        if profile_def{2} <= 0
            error('For {''REST'', duration_hours}, duration_hours must be positive.');
        end
        profile_name = 'REST';
        rest_hours_single = profile_def{2};
    else
        error(['sim_profile_LGM50 supports one dynamic profile token, ', ...
               '''CC'', or ''REST'' for single-profile simulation.']);
    end
end

function names = normalize_name_list(x)
    if isempty(x)
        names = {};
        return;
    end
    if ischar(x) || isstring(x)
        names = {char(x)};
        return;
    end
    if ~iscell(x)
        error('SPMeOutputs must be a string or cell array of strings.');
    end
    names = cell(size(x));
    for i = 1:numel(x)
        if ~(ischar(x{i}) || isstring(x{i}))
            error('SPMeOutputs entries must be strings.');
        end
        names{i} = char(x{i});
    end
end

function p = apply_spme_param_overrides(p, overrides)
    if isempty(overrides)
        return;
    end
    if ~isstruct(overrides)
        error('SPMeParamOverrides must be a struct when provided.');
    end
    names = fieldnames(overrides);
    for i = 1:numel(names)
        p.(names{i}) = overrides.(names{i});
    end
end

function diagnostics_file = save_spme_step_failure_diagnostics(diagnostics_file, context, ME, ...
        profile_name, init_soc, p, solver, t, t_profile, I, x, y, step_idx, ...
        req_spme_out, spme_out_vals)
    diagnostics_file = char_or_empty(diagnostics_file);
    if isempty(diagnostics_file)
        return;
    end

    diagnostic = struct();
    diagnostic.saved_at = datestr(now, 31);
    diagnostic.kind = 'spme_step_failure';
    diagnostic.context = context;
    diagnostic.error = exception_to_struct(ME);
    diagnostic.profile_name = profile_name;
    diagnostic.init_soc_input = init_soc;
    diagnostic.solver = solver;
    diagnostic.step_idx = step_idx;
    diagnostic.next_step_idx = step_idx + 1;
    diagnostic.step_time_s = value_or_nan(t, step_idx);
    diagnostic.next_step_time_s = value_or_nan(t, step_idx + 1);
    diagnostic.profile_time_s = value_or_nan(t_profile, step_idx);
    diagnostic.current_A = value_or_nan(I, step_idx);
    diagnostic.next_current_A = value_or_nan(I, step_idx + 1);
    diagnostic.previous_current_A = value_or_nan(I, step_idx - 1);
    diagnostic.p = p;
    diagnostic.state_indices = spme_state_indices(p);
    diagnostic.x0 = x(:, step_idx);
    diagnostic.y0 = y(:, step_idx);
    diagnostic.x0_summary = summarize_spme_state(diagnostic.x0, p);
    diagnostic.y0_summary = summarize_vector(diagnostic.y0);
    diagnostic.nonfinite_x0_idx = find(~isfinite(diagnostic.x0));
    diagnostic.nonfinite_y0_idx = find(~isfinite(diagnostic.y0));
    diagnostic.recent_window = recent_spme_window(t, I, x, y, step_idx, p, req_spme_out, spme_out_vals);
    diagnostic.partial_history = partial_spme_history(t, t_profile, I, x, y, ...
        step_idx, p, req_spme_out, spme_out_vals);

    [diagnostics_dir, ~, ~] = fileparts(diagnostics_file);
    if ~isempty(diagnostics_dir) && ~isfolder(diagnostics_dir)
        mkdir(diagnostics_dir);
    end
    save(diagnostics_file, 'diagnostic', '-v7.3');
end

function idx = spme_state_indices(p)
    idx = struct();
    idx.c_s_n = 1:(p.Nr - 1);
    idx.c_s_p = p.Nr:(2 * (p.Nr - 1));
    ce_start = 2 * p.Nr - 1;
    ce_end = ce_start + p.Nx - 4;
    idx.c_e = ce_start:ce_end;
    next_idx = ce_end + 1;
    if isfield(p, 'degradation_model') && ~strcmpi(char(p.degradation_model), 'none')
        idx.delta_sei = next_idx;
        next_idx = next_idx + 1;
    end
    if isfield(p, 'enable_lithium_plating') && logical(p.enable_lithium_plating)
        idx.delta_li = next_idx;
        next_idx = next_idx + 1;
    end
    if isfield(p, 'enable_anode_porosity_loss') && logical(p.enable_anode_porosity_loss)
        idx.epsilon_e_n = next_idx;
        next_idx = next_idx + 1;
    end
    if isfield(p, 'degradation_model') && strcmpi(char(p.degradation_model), 'safari')
        idx.c_ec_sei = next_idx:(next_idx + p.Nsei - 1);
    end
end

function summary = summarize_spme_state(xv, p)
    idx = spme_state_indices(p);
    summary = struct();
    summary.all = summarize_vector(xv);
    names = fieldnames(idx);
    for i = 1:numel(names)
        name = names{i};
        ii = idx.(name);
        ii = ii(ii >= 1 & ii <= numel(xv));
        if isempty(ii)
            summary.(name) = summarize_vector([]);
        else
            summary.(name) = summarize_vector(xv(ii));
        end
    end
end

function summary = summarize_vector(v)
    v = full(v(:));
    finite_v = v(isfinite(v));
    summary = struct();
    summary.numel = numel(v);
    summary.nonfinite_count = numel(v) - numel(finite_v);
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

function recent = recent_spme_window(t, I, x, y, step_idx, p, req_spme_out, spme_out_vals)
    first_idx = max(1, step_idx - 25);
    cols = first_idx:step_idx;
    recent = struct();
    recent.step_idx = cols(:);
    recent.time_s = t(cols);
    recent.current_A = I(cols);
    recent.y = y(:, cols);
    recent.x_summary = repmat(summarize_spme_state(x(:, cols(1)), p), numel(cols), 1);
    for k = 1:numel(cols)
        recent.x_summary(k) = summarize_spme_state(x(:, cols(k)), p);
    end
    if ~isempty(req_spme_out) && ~isempty(spme_out_vals)
        recent.spme_outputs = struct();
        for i = 1:numel(req_spme_out)
            recent.spme_outputs.(req_spme_out{i}) = spme_out_vals(i, cols).';
        end
    end
end

function hist = partial_spme_history(t, t_profile, I, x, y, step_idx, p, req_spme_out, spme_out_vals)
    cols = 1:step_idx;
    idx = spme_state_indices(p);
    hist = struct();
    hist.step_idx = cols(:);
    hist.time_s = t(cols);
    hist.profile_time_s = t_profile(cols);
    hist.current_A = I(cols);
    hist.voltage = y(1, cols).';
    if size(y, 1) >= 2
        hist.soc = y(2, cols).';
    else
        hist.soc = nan(numel(cols), 1);
    end
    hist.x = x(:, cols);
    hist.c_e_min = min(x(idx.c_e, cols), [], 1).';
    hist.c_e_mean = mean(x(idx.c_e, cols), 1).';
    hist.c_e_max = max(x(idx.c_e, cols), [], 1).';
    boundary = electrolyte_boundary_history(x(idx.c_e, cols), p);
    hist.ce0n = boundary.ce0n;
    hist.cens = boundary.cens;
    hist.cesp = boundary.cesp;
    hist.ce0p = boundary.ce0p;
    hist.electrolyte_log_arg_min = min([hist.ce0n, hist.cens, hist.cesp, hist.ce0p], [], 2);
    if isfield(idx, 'epsilon_e_n')
        hist.epsilon_e_n_state = x(idx.epsilon_e_n, cols).';
    else
        hist.epsilon_e_n_state = p.epsilon_e_n * ones(numel(cols), 1);
    end
    if isfield(idx, 'delta_sei')
        hist.delta_sei_state = x(idx.delta_sei, cols).';
    end
    if isfield(idx, 'delta_li')
        hist.delta_li_state = x(idx.delta_li, cols).';
    end
    if isfield(idx, 'c_ec_sei')
        hist.c_ec_sei_min = min(x(idx.c_ec_sei, cols), [], 1).';
        hist.c_ec_sei_mean = mean(x(idx.c_ec_sei, cols), 1).';
        hist.c_ec_sei_max = max(x(idx.c_ec_sei, cols), [], 1).';
    end
    if ~isempty(req_spme_out) && ~isempty(spme_out_vals)
        hist.spme_outputs = struct();
        for i = 1:numel(req_spme_out)
            hist.spme_outputs.(req_spme_out{i}) = spme_out_vals(i, cols).';
        end
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

function s = exception_to_struct(ME)
    s = struct();
    s.identifier = ME.identifier;
    s.message = ME.message;
    s.stack = ME.stack;
end

function v = value_or_nan(x, idx)
    if idx >= 1 && idx <= numel(x)
        v = x(idx);
    else
        v = NaN;
    end
end

function y = char_or_empty(x)
    if isempty(x)
        y = '';
    else
        y = char(x);
    end
end


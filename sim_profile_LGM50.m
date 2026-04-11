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
%                    {'none','ramadass','yang'} (default 'none')
%     'ElectrolyteModel' : SPMe electrolyte dynamics model
%                    {'fdm','fvm'} (default 'fdm')
%     'InitialState' : state struct for continuation from prior segment.
%                    For SPMe: struct('x', x_last)
%                    For DFN : struct('x', x_last, 'z', z_last)
%                    For ECM : struct('x', x_last)
%     'SuppressSummary' : true/false to suppress console summary print
%                    (default false)
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

    if ~isscalar(init_soc) || init_soc < 0 || init_soc > 1
        error('init_soc must be a scalar in [0, 1].');
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
    addParameter(parser, 'VEnd', []);
    addParameter(parser, 'InitialState', struct());
    addParameter(parser, 'SuppressSummary', false);
    addParameter(parser, 'SPMeOutputs', {});
    parse(parser, varargin{:});
    opt = parser.Results;

    base_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(base_dir, 'param'));
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'DFN'));
    addpath(fullfile(base_dir, 'ECM'));

    if ~isscalar(opt.Crate) || ~isfinite(opt.Crate)
        error('Crate must be a finite scalar.');
    end
    degradation_model = lower(char(opt.DegradationModel));
    if ~ismember(degradation_model, {'none', 'ramadass', 'yang'})
        error('DegradationModel must be one of: none, ramadass, yang.');
    end
    electrolyte_model = lower(char(opt.ElectrolyteModel));
    if ~ismember(electrolyte_model, {'fdm', 'fvm'})
        error('ElectrolyteModel must be one of: fdm, fvm.');
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
    calc = calculate_equilibrium_voltage_from_soc(init_soc);

    switch model_name
        case 'SPME'
            p = load_SPMe_params(init_soc, t_amb, Nr, delta_t, Nxn, Nxs, Nxp, disable_thermal, degradation_model, electrolyte_model);
            if is_cc_profile
                I = (p.OneC * opt.Crate) * ones(NT, 1);
            end

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
            nx_expected = numel(x0);
            if isstruct(opt.InitialState) && isfield(opt.InitialState, 'x') && ~isempty(opt.InitialState.x)
                x0 = opt.InitialState.x(:);
                if numel(x0) ~= nx_expected
                    error('InitialState.x size mismatch for SPMe continuation.');
                end
            end
            p.x0 = x0;

            spme_ss = SPMe(p, 'idas');
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
                [x(:, k+1), y(:, k+1)] = spme_ss.step(x(:, k), I(k));
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
            sim.c_e = x(2*p.Nr-1:2*p.Nr-1+p.Nx-4, 1:t_last);
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
            p = load_DFN_params(init_soc, t_amb, delta_t, disable_thermal, Nr, Nxn, Nxs, Nxp);
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

            p.x0 = [c_s_n; c_s_p; c_e; T];
            p.z0 = [Uref_n0; Uref_p0; i_en; i_ep; phi_e; jn; jp];
            [~, ~, p.y0] = dae_dfn(p.x0, p.z0, 0, p);
            if isstruct(opt.InitialState) && ~isempty(opt.InitialState)
                if isfield(opt.InitialState, 'x') && ~isempty(opt.InitialState.x)
                    p.x0 = opt.InitialState.x(:);
                    if numel(p.x0) ~= numel([c_s_n; c_s_p; c_e; T])
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

            [c_s_n_flat, c_s_p_flat, c_e_hist, T_hist] = parse_x(x, t_last, p);
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


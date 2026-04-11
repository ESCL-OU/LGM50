function results = test_ocp_vs_lli_3cycles(p, x0_fresh, lli_fracs, varargin)
%TEST_OCP_VS_LLI_3CYCLES
% Multi-cycle diagnostic for observing how imposed initial LLI propagates
% through the current SPMe implementation over 3 charge/discharge cycles.
%
% INPUTS
%   p          : (optional) parameter struct used by SPMe
%   x0_fresh   : (optional) valid SPMe state vector at rest / near-equilibrium
%   lli_fracs  : (optional) imposed initial LLI fractions, e.g. [0 0.02 0.05 0.1]
%
% OPTIONAL NAME-VALUE
%   'solver'        : 'rk' default
%   'C_rate'        : default 1/50
%   'Vmax'          : default from params (or 4.2 fallback)
%   'Vmin'          : default from params (or 2.8 fallback)
%   'rest_steps'    : default 0
%   'max_steps'     : default 200000
%   'n_cycles'      : default 3
%   'start_mode'    : 'charge' default
%   'useElectrolyteInventory' : true default
%   'init_soc'      : [] default; auto 0 for charge-start, 1 for discharge-start
%   'Tamb'          : default 25 [degC]
%   'DeltaT'        : default 1 [s]
%   'Nr'            : default 20
%   'Nxn'           : default 20
%   'Nxs'           : default 10
%   'Nxp'           : default 20
%   'DisableThermal': default true
%   'ElectrolyteModel' : default 'fvm'
%
% OUTPUT
%   results(k) for each imposed LLI case
%
% NOTES
%   - Uses the same OCP functions as the current SPMe code:
%       refPotentialAnode(...)
%       refPotentialCathode(...)
%   - This test is self-contained: if p/x0_fresh/lli_fracs are omitted,
%     they are generated internally from load_SPMe_params.
%   - Initial LLI is imposed diagnostically by removing Li uniformly from
%     the negative solid concentration profile only.
%   - Then the code is allowed to evolve normally over multiple cycles.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    ip = inputParser;
    addParameter(ip, 'solver', 'idas');
    addParameter(ip, 'C_rate', 1/2);
    addParameter(ip, 'Vmax', []);
    addParameter(ip, 'Vmin', []);
    addParameter(ip, 'rest_steps', 0);
    addParameter(ip, 'max_steps', 200000);
    addParameter(ip, 'n_cycles', 2);
    addParameter(ip, 'start_mode', 'charge');
    addParameter(ip, 'useElectrolyteInventory', true);
    addParameter(ip, 'ShowProgress', true);
    addParameter(ip, 'ProgressUpdateEvery', 2000);
    addParameter(ip, 'init_soc', 1.0);
    addParameter(ip, 'Tamb', 25);
    addParameter(ip, 'DeltaT', 1);
    addParameter(ip, 'Nr', 15);
    addParameter(ip, 'Nxn', 10);
    addParameter(ip, 'Nxs', 5);
    addParameter(ip, 'Nxp', 10);
    addParameter(ip, 'DisableThermal', true);
    addParameter(ip, 'ElectrolyteModel', 'fvm');
    parse(ip, varargin{:});
    opt = ip.Results;

    if nargin < 3 || isempty(lli_fracs)
        lli_fracs = [0 0.025 0.05 0.1];
    end
    lli_fracs = lli_fracs(:).';

    if nargin < 1 || isempty(p)
        init_soc_auto = opt.init_soc;
        p = load_SPMe_params(init_soc_auto, opt.Tamb, opt.Nr, opt.DeltaT, ...
            opt.Nxn, opt.Nxs, opt.Nxp, opt.DisableThermal, 'none', opt.ElectrolyteModel);
        fprintf('Initialized p internally via load_SPMe_params (init_soc=%.3f).\n', init_soc_auto);
    end
    if nargin < 2 || isempty(x0_fresh)
        if ~isfield(p, 'x0') || isempty(p.x0)
            error('x0_fresh was not provided and p.x0 is unavailable.');
        end
        x0_fresh = p.x0;
        fprintf('Using p.x0 as x0_fresh.\n');
    end

    if isempty(opt.Vmax)
        if isfield(p, 'volt_max')
            opt.Vmax = p.volt_max;
        else
            opt.Vmax = 4.2;
        end
    end
    if isempty(opt.Vmin)
        if isfield(p, 'volt_min')
            opt.Vmin = p.volt_min;
        else
            opt.Vmin = 2.8;
        end
    end

    % Run with degradation disabled during the diagnostic cycle sequence.
    % We want to observe propagation of an imposed initial inventory deficit.
    p_run = p;
    p_run.degradation_model = 'none';
    model = SPMe(p_run, opt.solver);
    fprintf('Using solver: %s\n', opt.solver);

    eps_s_n = get_eps_s_from_as(p_run.a_s_n, p_run.R_s_n);
    eps_s_p = get_eps_s_from_as(p_run.a_s_p, p_run.R_s_p);

    Q_nom_Ah = estimate_nominal_capacity_Ah(p_run, x0_fresh, eps_s_n, eps_s_p);
    I_slow = opt.C_rate * Q_nom_Ah;

    fprintf('Estimated nominal capacity = %.4f Ah\n', Q_nom_Ah);
    fprintf('Slow current = %.6f A (%.4f C)\n', I_slow, opt.C_rate);

    fresh = unpack_state(x0_fresh, p_run);
    nLi0 = compute_total_li_inventory(fresh, p_run, eps_s_n, eps_s_p, opt.useElectrolyteInventory);

    results = repmat(struct(), numel(lli_fracs), 1);
    total_legs = numel(lli_fracs) * opt.n_cycles * 2;
    completed_legs = 0;
    h_wait = [];
    if opt.ShowProgress
        try
            h_wait = waitbar(0, sprintf('Running OCP vs LLI test (0/%d legs)...', total_legs), ...
                'Name', 'test_ocp_vs_lli_3cycles');
            % Enlarge progress window for readability.
            try
                pos = get(h_wait, 'Position');
                pos(3) = max(pos(3), 700);
                pos(4) = max(pos(4), 140);
                set(h_wait, 'Position', pos);
            catch
                % no-op
            end
        catch
            h_wait = [];
        end
    end
    cleanup_wait = onCleanup(@() close_waitbar_safe(h_wait)); %#ok<NASGU>

    for k = 1:numel(lli_fracs)
        lli = lli_fracs(k);

        [x_init, lli_info] = impose_negative_solid_lli(x0_fresh, p_run, lli, eps_s_n, nLi0.total);

        x = x_init;
        for rr = 1:opt.rest_steps
            x = safe_integrate_state(model, x, 0);
        end

        legs = struct([]);
        global_idx = 0;
        t_global = 0;
        Ah_global = 0;

        start_with_charge = strcmpi(opt.start_mode, 'charge');

        for cyc = 1:opt.n_cycles
            if start_with_charge
                seq = {'charge','discharge'};
            else
                seq = {'discharge','charge'};
            end

            for jj = 1:2
                mode = seq{jj};

                if strcmpi(mode, 'charge')
                    % Sign convention in this repo: negative current = charge.
                    I = -abs(I_slow);
                    Vstop = opt.Vmax;
                else
                    I = +abs(I_slow);
                    Vstop = opt.Vmin;
                end

                leg = run_leg_multicycle( ...
                    model, x, I, Vstop, mode, opt.max_steps, ...
                    p_run, eps_s_n, eps_s_p, nLi0.total, ...
                    opt.useElectrolyteInventory, ...
                    cyc, jj, global_idx, t_global, Ah_global, ...
                    h_wait, completed_legs, total_legs, k, lli, ...
                    opt.ProgressUpdateEvery);

                legs = [legs; leg]; %#ok<AGROW>
                completed_legs = completed_legs + 1;
                update_waitbar_safe(h_wait, completed_legs / max(total_legs,1), ...
                    sprintf('Completed %d/%d legs (LLI %.1f%%, cycle %d, %s)', ...
                    completed_legs, total_legs, 100*lli, cyc, mode));
                if leg.hit_max_steps
                    warning(['Leg reached max_steps=%d without voltage stop (LLI %.1f%%, ', ...
                        'cycle %d, mode %s).'], opt.max_steps, 100*lli, cyc, mode);
                end

                x = leg.x_end;
                global_idx = leg.global_idx_end;
                t_global = leg.t_global_end;
                Ah_global = leg.Ah_global_end;

                for rr = 1:opt.rest_steps
                    x = safe_integrate_state(model, x, 0);
                    t_global = t_global + p_run.delta_t;
                end
            end
        end

        results(k).lli_frac = lli;
        results(k).lli_target_remove_mol = lli_info.n_remove_target_mol;
        results(k).lli_applied_remove_mol = lli_info.n_remove_applied_mol;
        results(k).lli_applied_frac_total = lli_info.n_remove_applied_mol / max(nLi0.total, eps);
        results(k).lli_reference_total_mol = nLi0.total;
        results(k).lli_negative_available_mol = lli_info.n_n_available_mol;
        results(k).lli_negative_available_frac_total = lli_info.n_n_available_mol / max(nLi0.total, eps);
        results(k).x_init = x_init;
        results(k).legs = legs;
        results(k).summary = summarize_case(legs);
    end

    fprintf('\n=== LLI request vs applied (negative-only imposition) ===\n');
    fprintf('Reference total Li inventory: %.6e mol\n', nLi0.total);
    for k = 1:numel(results)
        fprintf(['Case %d: requested=%.3f%% (%.6e mol), applied=%.3f%% (%.6e mol), ', ...
                 'negative-removable-limit=%.3f%% (%.6e mol)\n'], ...
            k, 100*results(k).lli_frac, results(k).lli_target_remove_mol, ...
            100*results(k).lli_applied_frac_total, results(k).lli_applied_remove_mol, ...
            100*results(k).lli_negative_available_frac_total, results(k).lli_negative_available_mol);
    end

    plot_results_multicycle(results, opt.n_cycles, opt.start_mode);
end


%% ========================= MULTICYCLE LEG =========================

function leg = run_leg_multicycle(model, x0, I, Vstop, mode, max_steps, ...
    p, eps_s_n, eps_s_p, nLi_ref_total, useElectrolyteInventory, ...
    cycle_number, leg_number, global_idx0, t_global0, Ah_global0, ...
    h_wait, completed_legs, total_legs, lli_case_idx, lli_frac, progress_update_every)

    x = x0;
    t_local = 0;
    t_global = t_global0;
    Ah_global = Ah_global0;
    dt = p.delta_t;

    arr_local_t  = nan(max_steps,1);
    arr_global_t = nan(max_steps,1);
    arr_V        = nan(max_steps,1);
    arr_Un       = nan(max_steps,1);
    arr_Up       = nan(max_steps,1);
    arr_Ucell    = nan(max_steps,1);
    arr_css_n    = nan(max_steps,1);
    arr_css_p    = nan(max_steps,1);
    arr_soc      = nan(max_steps,1);
    arr_Qleg_Ah  = nan(max_steps,1);
    arr_Qglob_Ah = nan(max_steps,1);
    arr_nLi_n    = nan(max_steps,1);
    arr_nLi_p    = nan(max_steps,1);
    arr_nLi_e    = nan(max_steps,1);
    arr_nLiTot   = nan(max_steps,1);
    arr_LLI      = nan(max_steps,1);
    arr_cycle    = nan(max_steps,1);
    arr_leg      = nan(max_steps,1);
    arr_I        = nan(max_steps,1);

    Qleg_Ah = 0;
    idx = 0;
    terminated_due_to_nan = false;
    hit_voltage_stop = false;
    termination_note = '';

    for k = 1:max_steps
        if mod(k, max(progress_update_every,1)) == 0
            frac = (completed_legs + min(k / max(max_steps,1), 0.999)) / max(total_legs, 1);
            update_waitbar_safe(h_wait, frac, sprintf(['Running leg %d/%d: case %d (LLI %.1f%%), ', ...
                'cycle %d, %s, step %d/%d'], ...
                completed_legs+1, total_legs, lli_case_idx, 100*lli_frac, ...
                cycle_number, mode, k, max_steps));
        end

        eval_on_next = false;
        try
            y = full(model.h('x', x, 'p', I).y);
            x_next = full(model.f_int('x0', x, 'p', I).x_next);
            x_eval = x;
        catch
            % Fallback for environments where direct access to CasADi
            % function-handle properties can fail indexing.
            [x_next, y] = model.step(x, I);
            eval_on_next = true;
            x_eval = x;
        end

        V = y(1);
        soc = y(2);

        st = unpack_state(x_eval, p);
        css_n = calc_css(st.c_s_n, I, p, 'n');
        css_p = calc_css(st.c_s_p, I, p, 'p');

        Un = refPotentialAnode(p, css_n / p.c_s_n_max);
        Up = refPotentialCathode(p, css_p / p.c_s_p_max);
        Ucell = Up - Un;

        nLi = compute_total_li_inventory(st, p, eps_s_n, eps_s_p, useElectrolyteInventory);

        if ~isfinite(V) || ~isfinite(soc) || ~isfinite(css_n) || ~isfinite(css_p)
            diag = state_numeric_diagnostic(st, p, I, css_n, css_p, V, soc);
            terminated_due_to_nan = true;
            termination_note = diag.message;
            warning('%s', diag.message);
            break;
        end

        idx = idx + 1;
        arr_local_t(idx)  = t_local;
        arr_global_t(idx) = t_global;
        arr_V(idx)        = V;
        arr_Un(idx)       = Un;
        arr_Up(idx)       = Up;
        arr_Ucell(idx)    = Ucell;
        arr_css_n(idx)    = css_n;
        arr_css_p(idx)    = css_p;
        arr_soc(idx)      = soc;
        arr_Qleg_Ah(idx)  = Qleg_Ah;
        arr_Qglob_Ah(idx) = Ah_global;
        arr_nLi_n(idx)    = nLi.n_n;
        arr_nLi_p(idx)    = nLi.n_p;
        arr_nLi_e(idx)    = nLi.n_e;
        arr_nLiTot(idx)   = nLi.total;
        arr_LLI(idx)      = 1 - nLi.total / nLi_ref_total;
        arr_cycle(idx)    = cycle_number;
        arr_leg(idx)      = leg_number;
        arr_I(idx)        = I;

        if strcmpi(mode, 'charge')
            if V >= Vstop
                hit_voltage_stop = true;
                if eval_on_next
                    x = x_next;
                end
                break;
            end
        else
            if V <= Vstop
                hit_voltage_stop = true;
                if eval_on_next
                    x = x_next;
                end
                break;
            end
        end

        x = x_next;

        t_local = t_local + dt;
        t_global = t_global + dt;
        dAh = abs(I) * dt / 3600;
        Qleg_Ah = Qleg_Ah + dAh;
        Ah_global = Ah_global + dAh;
    end

    leg.mode = mode;
    leg.cycle_number = cycle_number;
    leg.leg_number = leg_number;
    leg.x_end = x;

    leg.local_t_s = arr_local_t(1:idx);
    leg.global_t_s = arr_global_t(1:idx);
    leg.V = arr_V(1:idx);
    leg.Un = arr_Un(1:idx);
    leg.Up = arr_Up(1:idx);
    leg.Ucell = arr_Ucell(1:idx);
    leg.css_n = arr_css_n(1:idx);
    leg.css_p = arr_css_p(1:idx);
    leg.soc = arr_soc(1:idx);
    leg.Qleg_Ah = arr_Qleg_Ah(1:idx);
    leg.Qglob_Ah = arr_Qglob_Ah(1:idx);
    leg.nLi_n = arr_nLi_n(1:idx);
    leg.nLi_p = arr_nLi_p(1:idx);
    leg.nLi_e = arr_nLi_e(1:idx);
    leg.nLiTot = arr_nLiTot(1:idx);
    leg.LLI = arr_LLI(1:idx);
    leg.I = arr_I(1:idx);

    leg.global_idx_end = global_idx0 + idx;
    leg.t_global_end = t_global;
    leg.Ah_global_end = Ah_global;
    leg.hit_max_steps = (~terminated_due_to_nan) && (~hit_voltage_stop) && (idx >= max_steps);
    leg.terminated_due_to_nan = terminated_due_to_nan;
    leg.termination_note = termination_note;
end

function x_next = safe_integrate_state(model, x, I)
    try
        x_next = full(model.f_int('x0', x, 'p', I).x_next);
    catch
        [x_next, ~] = model.step(x, I);
    end
end

function update_waitbar_safe(h_wait, frac, msg)
    if isempty(h_wait) || ~ishandle(h_wait)
        return;
    end
    frac = min(max(frac, 0), 1);
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

function diag = state_numeric_diagnostic(st, p, I, css_n, css_p, V, soc)
    c_en = st.c_e(1:(p.Nxn-1));
    c_ep = st.c_e((p.Nxn-1)+p.Nxs : end);
    ce_n_edge = c_en(1);
    ce_p_edge = c_ep(end);
    theta_n = css_n / p.c_s_n_max;
    theta_p = css_p / p.c_s_p_max;
    i0_arg_n = (p.c_s_n_max - css_n) * css_n * ce_n_edge;
    i0_arg_p = (p.c_s_p_max - css_p) * css_p * ce_p_edge;
    msg = sprintf([ ...
        'Non-finite observable detected. I=%.6e A, V=%g, SoC=%g, ', ...
        'css_n=%g (theta_n=%g), css_p=%g (theta_p=%g), ', ...
        'ce_n_edge=%g, ce_p_edge=%g, i0_arg_n=%g, i0_arg_p=%g, ', ...
        'min(c_e)=%g, max(c_e)=%g.'], ...
        I, V, soc, css_n, theta_n, css_p, theta_p, ...
        ce_n_edge, ce_p_edge, i0_arg_n, i0_arg_p, ...
        min(st.c_e), max(st.c_e));
    diag = struct();
    diag.message = msg;
end


%% ========================= SUMMARY =========================

function S = summarize_case(legs)
    S = struct();

    nlegs = numel(legs);
    ncycles = nlegs / 2;

    S.charge_capacity_Ah = nan(ncycles,1);
    S.discharge_capacity_Ah = nan(ncycles,1);

    S.end_charge_nLi_n = nan(ncycles,1);
    S.end_charge_nLi_p = nan(ncycles,1);
    S.end_discharge_nLi_n = nan(ncycles,1);
    S.end_discharge_nLi_p = nan(ncycles,1);

    for c = 1:ncycles
        i1 = 2*c - 1;
        i2 = 2*c;

        if strcmpi(legs(i1).mode, 'charge')
            Lc = legs(i1);
            Ld = legs(i2);
        else
            Ld = legs(i1);
            Lc = legs(i2);
        end

        if ~isempty(Lc.Qleg_Ah)
            S.charge_capacity_Ah(c) = Lc.Qleg_Ah(end);
        else
            S.charge_capacity_Ah(c) = NaN;
        end
        if ~isempty(Ld.Qleg_Ah)
            S.discharge_capacity_Ah(c) = Ld.Qleg_Ah(end);
        else
            S.discharge_capacity_Ah(c) = NaN;
        end

        if ~isempty(Lc.nLi_n), S.end_charge_nLi_n(c) = Lc.nLi_n(end); else, S.end_charge_nLi_n(c) = NaN; end
        if ~isempty(Lc.nLi_p), S.end_charge_nLi_p(c) = Lc.nLi_p(end); else, S.end_charge_nLi_p(c) = NaN; end
        if ~isempty(Ld.nLi_n), S.end_discharge_nLi_n(c) = Ld.nLi_n(end); else, S.end_discharge_nLi_n(c) = NaN; end
        if ~isempty(Ld.nLi_p), S.end_discharge_nLi_p(c) = Ld.nLi_p(end); else, S.end_discharge_nLi_p(c) = NaN; end
    end
end


%% ========================= STATE / INVENTORY HELPERS =========================

function st = unpack_state(x, p)
    idx = 0;

    st.c_s_n = x(idx + (1:(p.Nr-1))); idx = idx + (p.Nr-1);
    st.c_s_p = x(idx + (1:(p.Nr-1))); idx = idx + (p.Nr-1);
    st.c_e   = x(idx + (1:(p.Nx-3))); idx = idx + (p.Nx-3);

    if length(x) > idx
        st.delta_sei = x(idx+1);
    else
        st.delta_sei = [];
    end
end

function x = pack_state(st, p)
    if isfield(st, 'delta_sei') && ~isempty(st.delta_sei)
        x = [st.c_s_n; st.c_s_p; st.c_e; st.delta_sei];
    else
        x = [st.c_s_n; st.c_s_p; st.c_e];
    end
end

function eps_s = get_eps_s_from_as(a_s, R_s)
    eps_s = a_s * R_s / 3;
end

function A = get_cell_area(p)
    if isfield(p, 'Area')
        A = p.Area;
    elseif isfield(p, 'Area_n')
        A = p.Area_n;
    elseif isfield(p, 'Area_p')
        A = p.Area_p;
    else
        error('Could not find p.Area or p.Area_n / p.Area_p');
    end
end

function css = calc_css(c_s_interior, I, p, whichEl)
    [~, ~, ~, ~, C_n, C_p, D_n, D_p] = spm_plant_obs_mats(p);

    if strcmpi(whichEl, 'n')
        jn = I / (p.Faraday * p.a_s_n * p.Area_n * p.L_n);
        css = full(C_n * c_s_interior + D_n * jn);
    else
        jp = -I / (p.Faraday * p.a_s_p * p.Area_p * p.L_p);
        css = full(C_p * c_s_interior + D_p * jp);
    end
end

function c_full = reconstruct_full_radial_profile(c_interior, I, p, whichEl)
    [~, ~, ~, ~, C_n, C_p, D_n, D_p] = spm_plant_obs_mats(p);

    if strcmpi(whichEl, 'n')
        j = I / (p.Faraday * p.a_s_n * p.Area_n * p.L_n);
        css = full(C_n * c_interior + D_n * j);
    else
        j = -I / (p.Faraday * p.a_s_p * p.Area_p * p.L_p);
        css = full(C_p * c_interior + D_p * j);
    end

    c0 = c_interior(1);
    c_full = [c0; c_interior; css];
end

function cavg = radial_average(c_full)
    Nr_full = length(c_full);
    r = linspace(0, 1, Nr_full).';
    num = trapz(r, c_full .* (r.^2));
    den = trapz(r, r.^2);
    cavg = num / den;
end

function nLi = compute_total_li_inventory(st, p, eps_s_n, eps_s_p, includeElectrolyte)
    c_n_full = reconstruct_full_radial_profile(st.c_s_n, 0, p, 'n');
    c_p_full = reconstruct_full_radial_profile(st.c_s_p, 0, p, 'p');

    cavg_n = radial_average(c_n_full);
    cavg_p = radial_average(c_p_full);

    Vsn = eps_s_n * p.Area_n * p.L_n;
    Vsp = eps_s_p * p.Area_p * p.L_p;

    n_n = cavg_n * Vsn;
    n_p = cavg_p * Vsp;

    n_e = 0;
    if includeElectrolyte
        % Use interior electrolyte states directly to avoid dependency on
        % boundary reconstruction helper signatures.
        c_en = st.c_e(1:(p.Nxn-1));
        c_es = st.c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
        c_ep = st.c_e((p.Nxn-1)+p.Nxs : end);

        ce_n_bar = mean(c_en);
        ce_s_bar = mean(c_es);
        ce_p_bar = mean(c_ep);

        n_e = ...
            ce_n_bar * p.epsilon_e_n * p.Area_n * p.L_n + ...
            ce_s_bar * p.epsilon_e_s * p.Area_s * p.L_s + ...
            ce_p_bar * p.epsilon_e_p * p.Area_p * p.L_p;
    end

    nLi.n_n = n_n;
    nLi.n_p = n_p;
    nLi.n_e = n_e;
    nLi.total = n_n + n_p + n_e;
end

function [x_new, info] = impose_negative_solid_lli(x_old, p, lli_frac, eps_s_n, nLi_ref_total)
    st = unpack_state(x_old, p);

    c_n_full = reconstruct_full_radial_profile(st.c_s_n, 0, p, 'n');
    Vsn = eps_s_n * p.Area_n * p.L_n;

    % IMPORTANT: lli_frac is defined against TOTAL cell Li inventory.
    n_remove_target = lli_frac * nLi_ref_total;
    n_remove = n_remove_target;
    % No capping: apply the requested LLI amount directly to negative solid.
    dc_n = n_remove / max(Vsn, eps);

    c_n_full_new = c_n_full - dc_n;

    st.c_s_n = c_n_full_new(2:end-1);
    x_new = pack_state(st, p);

    info = struct();
    info.n_remove_target_mol = n_remove_target;
    info.n_remove_applied_mol = n_remove;
    info.n_n_available_mol = NaN;
end

function Qnom_Ah = estimate_nominal_capacity_Ah(p, x0, eps_s_n, eps_s_p)
    % Prefer model-calibrated 1C definition when available.
    if isfield(p, 'OneC') && isscalar(p.OneC) && isfinite(p.OneC) && p.OneC > 0
        Qnom_Ah = p.OneC;
        return;
    end

    st = unpack_state(x0, p);

    c_n_full = reconstruct_full_radial_profile(st.c_s_n, 0, p, 'n');
    c_p_full = reconstruct_full_radial_profile(st.c_s_p, 0, p, 'p');

    cavg_n = radial_average(c_n_full);
    cavg_p = radial_average(c_p_full);

    Vsn = eps_s_n * p.Area_n * p.L_n;
    Qnom_Ah = p.Faraday * Vsn * max(1e-6, min(cavg_n, p.c_s_n_max - cavg_n)) / 3600;

    if Qnom_Ah < 0.01
        Vsp = eps_s_p * p.Area_p * p.L_p;
        Qnom_Ah = p.Faraday * Vsp * max(1e-6, min(cavg_p, p.c_s_p_max - cavg_p)) / 3600;
    end

    Qnom_Ah = max(Qnom_Ah, 0.1);
end


%% ========================= PLOTTING =========================

function plot_results_multicycle(results, n_cycles, start_mode)
    ncase = numel(results);
    styles = make_case_styles(ncase);
    % Keep cycle-overlay layout fixed for readability:
    % left column = discharge, right column = charge.
    mode_a = 'discharge';
    mode_b = 'charge';

    figure('Name','Global evolution over 3 cycles','Color','w');
    tiledlayout(3,2);

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, Uc] = concat_field(results(k).legs, 'Qglob_Ah', 'Ucell');
        plot_case(Qg, Uc, styles(k), case_label(results(k)));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Cell OCP U_p-U_n [V]');
    title('Global cell OCP');
    legend('Location','best');

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, V] = concat_field(results(k).legs, 'Qglob_Ah', 'V');
        plot_case(Qg, V, styles(k), case_label(results(k)));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Terminal voltage [V]');
    title('Global terminal voltage');
    legend('Location','best');

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, Un] = concat_field(results(k).legs, 'Qglob_Ah', 'Un');
        plot_case(Qg, Un, styles(k), case_label(results(k)));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Anode OCP [V]');
    title('Global anode OCP');
    legend('Location','best');

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, Up] = concat_field(results(k).legs, 'Qglob_Ah', 'Up');
        plot_case(Qg, Up, styles(k), case_label(results(k)));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Cathode OCP [V]');
    title('Global cathode OCP');
    legend('Location','best');

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, nN] = concat_field(results(k).legs, 'Qglob_Ah', 'nLi_n');
        plot_case(Qg, nN, styles(k), sprintf('Neg %s', case_label(results(k))));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Negative solid Li [mol]');
    title('Negative inventory');
    legend('Location','best');

    nexttile; hold on; grid on;
    for k = 1:ncase
        [Qg, nP] = concat_field(results(k).legs, 'Qglob_Ah', 'nLi_p');
        plot_case(Qg, nP, styles(k), sprintf('Pos %s', case_label(results(k))));
    end
    xlabel('Cumulative throughput [Ah]');
    ylabel('Positive solid Li [mol]');
    title('Positive inventory');
    legend('Location','best');


    % Per-cycle overlays
    for cyc = 1:n_cycles
        figure('Name', sprintf('Cycle %d overlays', cyc), 'Color', 'w');
        tiledlayout(3,2);

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, Uc] = get_cycle_leg_overlay(results(k).legs, cyc, mode_a, 'Qleg_Ah', 'Ucell');
            plot_case(Q, Uc, styles(k), case_label(results(k)));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_a));
        ylabel('Cell OCP [V]');
        title(sprintf('Cycle %d %s: cell OCP', cyc, mode_a));
        legend('Location','best');

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, Uc] = get_cycle_leg_overlay(results(k).legs, cyc, mode_b, 'Qleg_Ah', 'Ucell');
            plot_case(Q, Uc, styles(k), case_label(results(k)));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_b));
        ylabel('Cell OCP [V]');
        title(sprintf('Cycle %d %s: cell OCP', cyc, mode_b));
        legend('Location','best');

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, Up] = get_cycle_leg_overlay(results(k).legs, cyc, mode_a, 'Qleg_Ah', 'Up');
            st_pos = styles(k);
            st_neg = styles(k);
            st_neg.LineStyle = ':';
            st_neg.Marker = 'none';
            st_neg.MarkerIndices = [];
            [~, Un] = get_cycle_leg_overlay(results(k).legs, cyc, mode_a, 'Qleg_Ah', 'Un');
            plot_case(Q, Up, st_pos, sprintf('%s (U_p)', case_label(results(k))));
            plot_case(Q, Un, st_neg, sprintf('%s (U_n)', case_label(results(k))));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_a));
        ylabel('Electrode potential [V]');
        title(sprintf('Cycle %d %s: electrode potentials', cyc, mode_a));
        legend('Location','best');

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, Up] = get_cycle_leg_overlay(results(k).legs, cyc, mode_b, 'Qleg_Ah', 'Up');
            st_pos = styles(k);
            st_neg = styles(k);
            st_neg.LineStyle = ':';
            st_neg.Marker = 'none';
            st_neg.MarkerIndices = [];
            [~, Un] = get_cycle_leg_overlay(results(k).legs, cyc, mode_b, 'Qleg_Ah', 'Un');
            plot_case(Q, Up, st_pos, sprintf('%s (U_p)', case_label(results(k))));
            plot_case(Q, Un, st_neg, sprintf('%s (U_n)', case_label(results(k))));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_b));
        ylabel('Electrode potential [V]');
        title(sprintf('Cycle %d %s: electrode potentials', cyc, mode_b));
        legend('Location','best');

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, nP] = get_cycle_leg_overlay(results(k).legs, cyc, mode_a, 'Qleg_Ah', 'nLi_p');
            st_pos = styles(k);
            st_neg = styles(k);
            st_neg.LineStyle = ':';
            st_neg.Marker = 'none';
            st_neg.MarkerIndices = [];
            [~, nN] = get_cycle_leg_overlay(results(k).legs, cyc, mode_a, 'Qleg_Ah', 'nLi_n');
            plot_case(Q, nP, st_pos, sprintf('%s (pos)', case_label(results(k))));
            plot_case(Q, nN, st_neg, sprintf('%s (neg)', case_label(results(k))));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_a));
        ylabel('Solid Li inventory [mol]');
        title(sprintf('Cycle %d %s: pos + neg inventory', cyc, mode_a));
        legend('Location','best');

        nexttile; hold on; grid on;
        for k = 1:ncase
            [Q, nP] = get_cycle_leg_overlay(results(k).legs, cyc, mode_b, 'Qleg_Ah', 'nLi_p');
            st_pos = styles(k);
            st_neg = styles(k);
            st_neg.LineStyle = ':';
            st_neg.Marker = 'none';
            st_neg.MarkerIndices = [];
            [~, nN] = get_cycle_leg_overlay(results(k).legs, cyc, mode_b, 'Qleg_Ah', 'nLi_n');
            plot_case(Q, nP, st_pos, sprintf('%s (pos)', case_label(results(k))));
            plot_case(Q, nN, st_neg, sprintf('%s (neg)', case_label(results(k))));
        end
        xlabel(sprintf('%s capacity in leg [Ah]', mode_b));
        ylabel('Solid Li inventory [mol]');
        title(sprintf('Cycle %d %s: pos + neg inventory', cyc, mode_b));
        legend('Location','best');
    end

    % Summary printout
    fprintf('\n=== Capacity summary over cycles ===\n');
    for k = 1:ncase
        fprintf('LLI %.2f%%\n', 100*results(k).lli_frac);
        for c = 1:numel(results(k).summary.charge_capacity_Ah)
            fprintf('  Cycle %d: charge %.4f Ah, discharge %.4f Ah\n', ...
                c, ...
                results(k).summary.charge_capacity_Ah(c), ...
                results(k).summary.discharge_capacity_Ah(c));
        end
    end
end

function lbl = case_label(res_k)
    lbl = sprintf('LLI req %.1f%% / app %.2f%%', ...
        100*res_k.lli_frac, 100*res_k.lli_applied_frac_total);
end

function styles = make_case_styles(ncase)
    base_colors = lines(max(ncase, 7));
    line_styles = {'-', '--', ':', '-.'};
    markers = {'none', 'o', 's', 'd', '^', 'v', 'x', '+'};
    styles = repmat(struct('Color', [0 0 0], 'LineStyle', '-', ...
        'LineWidth', 1.4, 'Marker', 'none', 'MarkerIndices', []), ncase, 1);
    for k = 1:ncase
        styles(k).Color = base_colors(mod(k-1, size(base_colors,1)) + 1, :);
        styles(k).LineStyle = line_styles{mod(k-1, numel(line_styles)) + 1};
        styles(k).Marker = markers{mod(k-1, numel(markers)) + 1};
        if strcmp(styles(k).Marker, 'none')
            styles(k).MarkerIndices = [];
        else
            styles(k).MarkerIndices = 1:80:200000;
        end
    end
end

function plot_case(x, y, st, disp_name)
    if isempty(x) || isempty(y)
        return;
    end
    mask = isfinite(x) & isfinite(y);
    x = x(mask);
    y = y(mask);
    if isempty(x) || isempty(y)
        return;
    end
    % If only a few samples exist, force markers so the trace is visible.
    force_marker = numel(x) <= 2;
    if isempty(st.MarkerIndices) && ~force_marker
        plot(x, y, ...
            'Color', st.Color, ...
            'LineStyle', st.LineStyle, ...
            'LineWidth', st.LineWidth, ...
            'DisplayName', disp_name);
    else
        plot(x, y, ...
            'Color', st.Color, ...
            'LineStyle', st.LineStyle, ...
            'LineWidth', st.LineWidth, ...
            'Marker', choose_marker(st.Marker, force_marker), ...
            'MarkerIndices', choose_marker_indices(st.MarkerIndices, numel(x), force_marker), ...
            'MarkerSize', 6, ...
            'DisplayName', disp_name);
    end
end

function [xcat, ycat] = concat_field(legs, xfield, yfield)
    xcat = [];
    ycat = [];
    for i = 1:numel(legs)
        xcat = [xcat; legs(i).(xfield)(:)];
        ycat = [ycat; legs(i).(yfield)(:)];
    end
end

function [x, y] = get_cycle_leg_overlay(legs, cyc, mode, xfield, yfield)
    x = [];
    y = [];
    best_n = -1;
    best_x = [];
    best_y = [];
    for i = 1:numel(legs)
        if legs(i).cycle_number == cyc && strcmpi(legs(i).mode, mode)
            xi = legs(i).(xfield)(:);
            yi = legs(i).(yfield)(:);
            if isempty(xi) || isempty(yi)
                continue;
            end
            n = min(numel(xi), numel(yi));
            xi = xi(1:n);
            yi = yi(1:n);
            mask = isfinite(xi) & isfinite(yi);
            n_valid = nnz(mask);
            if n_valid > best_n
                best_n = n_valid;
                best_x = xi(mask);
                best_y = yi(mask);
            end
        end
    end
    x = best_x;
    y = best_y;
end

function m = choose_marker(marker_in, force_marker)
    if force_marker
        m = 'o';
    else
        m = marker_in;
    end
end

function idx = choose_marker_indices(marker_idx_in, npts, force_marker)
    if force_marker
        idx = 1:max(npts,1);
    else
        idx = marker_idx_in;
    end
end
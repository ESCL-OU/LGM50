function out = test_spme_electrolyte_degradation_mode_compare(profile_sequence, init_soc, electrolyte_model)
%TEST_SPME_ELECTROLYTE_DEGRADATION_MODE_COMPARE
% Compare electrolyte depletion under SEI/plating configurations.
%
% Usage:
%   out = test_spme_electrolyte_degradation_mode_compare();
%   out = test_spme_electrolyte_degradation_mode_compare({'CC', -1.0, 'CC', 1.0}, 0.1);
%   out = test_spme_electrolyte_degradation_mode_compare([], [], 'fvm');

    if nargin < 1 || isempty(profile_sequence)
        profile_sequence = {'REST', 0.25, 'CC', -0.2, 'REST', 0.25, 'CC', 0.2};
    end
    if nargin < 2 || isempty(init_soc)
        init_soc = 0.1;
    end
    if nargin < 3 || isempty(electrolyte_model)
        electrolyte_model = 'fdm';
    end
    electrolyte_model = lower(char(electrolyte_model));
    if ~ismember(electrolyte_model, {'fdm', 'fvm'})
        error('electrolyte_model must be ''fdm'' or ''fvm''.');
    end

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    delta_t = 1;
    t_amb = 25;
    Nr = 15;
    Nxn = 20;
    Nxs = 5;
    Nxp = 20;
    degradation_model = 'safari';
    repeat_count = 2;

    req_outputs = {'delta_sei', 'delta_sei_dot', ...
        'delta_li', 'delta_li_dot', 'delta_film', 'delta_film_dot', ...
        'epsilon_e_n', 'epsilon_e_n_dot', ...
        'j_sei', 'j_plating', 'j_side', 'R_film', ...
        'V_electrolytePolar', 'ce0n', 'cens', 'cesp', 'ce0p', ...
        'cen_bar', 'ces_bar', 'cep_bar'};

    cases = struct('label', {}, 'degradation_model', {}, ...
        'enable_lithium_plating', {}, 'enable_anode_porosity_loss', {}, ...
        'overrides', {});
    cases(end+1) = struct('label', 'No SEI / no plating', ...
        'degradation_model', 'none', ...
        'enable_lithium_plating', false, ...
        'enable_anode_porosity_loss', false, ...
        'overrides', struct());
    cases(end+1) = struct('label', 'SEI only (Safari)', ...
        'degradation_model', degradation_model, ...
        'enable_lithium_plating', false, ...
        'enable_anode_porosity_loss', true, ...
        'overrides', struct());
    cases(end+1) = struct('label', 'SEI only (Ramadass)', ...
        'degradation_model', 'ramadass', ...
        'enable_lithium_plating', false, ...
        'enable_anode_porosity_loss', true, ...
        'overrides', struct());
    cases(end+1) = struct('label', 'Li plating only (SEI k0=0)', ...
        'degradation_model', degradation_model, ...
        'enable_lithium_plating', true, ...
        'enable_anode_porosity_loss', true, ...
        'overrides', struct('k0_sei', 0, 'i0s', 0));
    cases(end+1) = struct('label', 'SEI + Li plating', ...
        'degradation_model', degradation_model, ...
        'enable_lithium_plating', true, ...
        'enable_anode_porosity_loss', true, ...
        'overrides', struct());

    traces = struct('label', {}, 'sim', {}, 'time_hr', {}, 'current', {}, ...
        'voltage', {}, 'soc', {}, 'c_e_min', {}, 'c_e_mean', {}, 'c_e_max', {}, ...
        'ce0n', {}, 'cens', {}, 'cesp', {}, 'ce0p', {}, 'log_arg_min', {}, ...
        'epsilon_e_n', {}, 'delta_sei', {}, 'delta_li', {}, 'delta_film', {}, ...
        'R_film', {}, 'j_sei', {}, 'j_plating', {}, ...
        'electrolyte_inventory', {}, 'electrolyte_inventory_rel_change', {});

    fprintf('Electrolyte comparison uses repeat_count=%d from fresh init_soc=%.12g.\n', ...
        repeat_count, init_soc);
    if strcmp(electrolyte_model, 'fdm')
        fprintf(['NOTE: FDM keeps the legacy non-conservative diffusion operator; ', ...
            'use electrolyte_model=''fvm'' for the conservation diagnostic.\n']);
    end
    for i = 1:numel(cases)
        cfg = cases(i);
        fprintf('Running case %d/%d: %s\n', i, numel(cases), cfg.label);
        sim = sim_profile_sequence_LGM50('SPMe', init_soc, profile_sequence, ...
            'Repeat', repeat_count, ...
            'DeltaT', delta_t, ...
            'Tamb', t_amb, ...
            'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
            'DegradationModel', cfg.degradation_model, ...
            'EnableLithiumPlating', cfg.enable_lithium_plating, ...
            'EnableAnodePorosityLoss', cfg.enable_anode_porosity_loss, ...
            'ElectrolyteModel', electrolyte_model, ...
            'SPMeOutputs', req_outputs, ...
            'SPMeParamOverrides', cfg.overrides, ...
            'MemoryLight', false, ...
            'ShowProgress', false, ...
            'SaveOutput', false);
        tr = trace_from_sim(sim, cfg.label);
        traces(end+1) = tr; %#ok<AGROW>
    end

    plot_comparison(traces);

    out = struct();
    out.profile_sequence = profile_sequence;
    out.init_soc = init_soc;
    out.electrolyte_model = electrolyte_model;
    out.note = ['FVM is recommended for electrolyte conservation checks. ', ...
        'FDM still uses the legacy non-conservative diffusion operator.'];
    out.repeat_count = repeat_count;
    out.cases = cases;
    out.traces = traces;
end

function tr = trace_from_sim(sim, label)
    outputs = get_outputs(sim);
    p = sim.p;
    x = sim.x;
    idx_ce_start = 2 * p.Nr - 1;
    idx_ce_end = idx_ce_start + p.Nx - 4;
    c_e = x(idx_ce_start:idx_ce_end, :);

    ce0n = get_output_or_nan(outputs, 'ce0n', numel(sim.time));
    cens = get_output_or_nan(outputs, 'cens', numel(sim.time));
    cesp = get_output_or_nan(outputs, 'cesp', numel(sim.time));
    ce0p = get_output_or_nan(outputs, 'ce0p', numel(sim.time));

    tr = struct();
    tr.label = label;
    tr.sim = sim;
    tr.time_hr = sim.time(:) / 3600;
    tr.current = sim.current(:);
    tr.voltage = sim.voltage(:);
    tr.soc = sim.soc(:);
    tr.c_e_min = min(c_e, [], 1).';
    tr.c_e_mean = mean(c_e, 1).';
    tr.c_e_max = max(c_e, [], 1).';
    tr.ce0n = ce0n;
    tr.cens = cens;
    tr.cesp = cesp;
    tr.ce0p = ce0p;
    tr.log_arg_min = min([ce0n, cens, cesp, ce0p], [], 2);
    tr.epsilon_e_n = get_output_or_nan(outputs, 'epsilon_e_n', numel(sim.time));
    if all(isnan(tr.epsilon_e_n))
        tr.epsilon_e_n = p.epsilon_e_n * ones(numel(sim.time), 1);
    end
    tr.delta_sei = get_output_or_nan(outputs, 'delta_sei', numel(sim.time));
    tr.delta_li = get_output_or_nan(outputs, 'delta_li', numel(sim.time));
    tr.delta_film = get_output_or_nan(outputs, 'delta_film', numel(sim.time));
    tr.R_film = get_output_or_nan(outputs, 'R_film', numel(sim.time));
    tr.j_sei = get_output_or_nan(outputs, 'j_sei', numel(sim.time));
    tr.j_plating = get_output_or_nan(outputs, 'j_plating', numel(sim.time));
    tr.electrolyte_inventory = electrolyte_inventory(c_e, tr.epsilon_e_n, p);
    tr.electrolyte_inventory_rel_change = ...
        (tr.electrolyte_inventory - tr.electrolyte_inventory(1)) / ...
        max(abs(tr.electrolyte_inventory(1)), eps);
end

function outputs = get_outputs(sim)
    if isfield(sim, 'spme_outputs')
        outputs = sim.spme_outputs;
    else
        outputs = struct();
    end
end

function y = get_output_or_nan(outputs, name, n)
    if isfield(outputs, name)
        y = outputs.(name)(:);
    else
        y = nan(n, 1);
    end
end

function plot_comparison(traces)
    colors = lines(numel(traces));

    figure('Name', 'SPMe electrolyte comparison: degradation modes', 'Color', 'w');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_field(traces, colors, 'current', 'Current [A]', 'Current');
    plot_field(traces, colors, 'voltage', 'Voltage [V]', 'Voltage');
    plot_field(traces, colors, 'c_e_min', 'min c_e [mol/m^3]', 'Internal Electrolyte Minimum');
    plot_field(traces, colors, 'c_e_mean', 'mean c_e [mol/m^3]', 'Internal Electrolyte Mean');
    plot_field(traces, colors, 'ce0n', 'ce0n [mol/m^3]', 'Anode Boundary Electrolyte');
    plot_field(traces, colors, 'log_arg_min', 'min log arg [mol/m^3]', 'Minimum Log Argument');

    figure('Name', 'SPMe electrolyte inventory comparison: degradation modes', 'Color', 'w');
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_field(traces, colors, 'electrolyte_inventory', 'N_e [mol]', ...
        'Weighted Electrolyte Lithium Inventory');
    plot_field(traces, colors, 'electrolyte_inventory_rel_change', '\DeltaN_e/N_e(0) [-]', ...
        'Relative Inventory Change');

    figure('Name', 'SPMe degradation comparison: side reactions', 'Color', 'w');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_field(traces, colors, 'epsilon_e_n', '\epsilon_{e,n} [-]', 'Anode Porosity');
    plot_field(traces, colors, 'delta_sei', '\delta_{SEI} [m]', 'SEI Thickness');
    plot_field(traces, colors, 'delta_li', '\delta_{Li} [m]', 'Lithium Plating Thickness');
    plot_field(traces, colors, 'delta_film', '\delta_{film} [m]', 'Total Film Thickness');
    plot_field(traces, colors, 'j_sei', 'j_{SEI} [mol/(m^2 s)]', 'SEI Flux');
    plot_field(traces, colors, 'j_plating', 'j_{plating} [mol/(m^2 s)]', 'Plating Flux');
end

function plot_field(traces, colors, field_name, y_label, plot_title)
    nexttile;
    hold on;
    for i = 1:numel(traces)
        plot(traces(i).time_hr, traces(i).(field_name), ...
            'Color', colors(i, :), 'LineWidth', 1.1);
    end
    if contains(y_label, 'mol/m^3')
        yline(0, ':r');
    end
    grid on;
    xlabel('Time [h]');
    ylabel(y_label);
    title(plot_title);
    legend({traces.label}, 'Location', 'best', 'Interpreter', 'none');
end

function inventory = electrolyte_inventory(c_e, epsilon_e_n, p)
    n_n = p.Nxn - 1;
    n_s = p.Nxs - 1;
    n_p = p.Nxp - 1;
    dx_n = p.L_n * p.delta_x_n;
    dx_s = p.L_s * p.delta_x_s;
    dx_p = p.L_p * p.delta_x_p;
    inventory = zeros(size(c_e, 2), 1);
    for k = 1:size(c_e, 2)
        eps_n = epsilon_e_n(k);
        if ~isfinite(eps_n)
            eps_n = p.epsilon_e_n;
        end
        weights = [eps_n * p.Area_n * dx_n * ones(n_n, 1); ...
                   p.epsilon_e_s * p.Area_s * dx_s * ones(n_s, 1); ...
                   p.epsilon_e_p * p.Area_p * dx_p * ones(n_p, 1)];
        inventory(k) = weights.' * c_e(:, k);
    end
end

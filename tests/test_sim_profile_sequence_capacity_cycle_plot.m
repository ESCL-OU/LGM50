function out = test_sim_profile_sequence_capacity_cycle_plot()
%TEST_SIM_PROFILE_SEQUENCE_CAPACITY_CYCLE_PLOT
% Validate capacity-cycle bookkeeping in sim_profile_sequence_LGM50 using:
%   {1/10C discharge, 1C charge, UDDS, 1/10C charge}
% starting from SoC = 100%, and plot charged/discharged Ah per cycle.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    model_name = 'SPMe';
    init_soc = 0.0;
    sequence = {'REST', 0.5, 'CC', -1.0, 'REST', 0.5, 'CC', 1.0, 'REST', 0.5, 'CC', -1.0, 'REST', 0.5, 'CC', 1.0};
    repeat_count = 2;
    degradation_model = 'none';

    sim = sim_profile_sequence_LGM50(model_name, init_soc, sequence, ...
        'Repeat', repeat_count, ...
        'DeltaT', 1, ...
        'Tamb', 25, ...
        'SPMeOutputs', {'delta_sei', 'cen_bar', 'ces_bar', 'cep_bar', ...
            'cs_ave_n', 'cs_ave_p', ...
            'c_en_dot_bar', 'c_es_dot_bar', 'c_ep_dot_bar', ...
            'cs_ave_n_dot', 'cs_ave_p_dot'}, ...
        'DegradationModel', degradation_model, ...
        'ElectrolyteModel', 'fvm', ...
        'SaveOutput', false, ...
        'ShowProgress', false);

    assert(isfield(sim, 'capacity_cycle_log'), ...
        'Missing sim.capacity_cycle_log.');
    assert(isfield(sim, 'capacity_cycle'), ...
        'Missing sim.capacity_cycle time-series label.');
    assert(isfield(sim, 'segment_log'), ...
        'Missing sim.segment_log.');

    cycle_log = sim.capacity_cycle_log;
    n_capacity_cycles = numel(cycle_log);
    expected_cycles = repeat_count * 2;
    assert(n_capacity_cycles == expected_cycles, ...
        'Unexpected number of capacity cycles. Expected %d, got %d.', ...
        expected_cycles, n_capacity_cycles);

    charged_ah = zeros(n_capacity_cycles, 1);
    discharged_ah = zeros(n_capacity_cycles, 1);
    net_ah = zeros(n_capacity_cycles, 1);
    cycle_id = zeros(n_capacity_cycles, 1);
    for k = 1:n_capacity_cycles
        cycle_id(k) = cycle_log(k).cycle;
        charged_ah(k) = cycle_log(k).charged_Ah;
        discharged_ah(k) = cycle_log(k).discharged_Ah;
        net_ah(k) = cycle_log(k).net_Ah;
        assert(charged_ah(k) >= 0 && discharged_ah(k) >= 0, ...
            'Cycle %d has negative charged/discharged capacity.', k);
    end

    cycle_edges = find(diff(sim.capacity_cycle) > 0) + 1;
    cycle_switch_t = sim.time(cycle_edges);

    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'delta_sei')
        delta_sei = sim.spme_outputs.delta_sei(:);
    else
        delta_sei = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cen_bar')
        cen_bar = sim.spme_outputs.cen_bar(:);
    else
        cen_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'ces_bar')
        ces_bar = sim.spme_outputs.ces_bar(:);
    else
        ces_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cep_bar')
        cep_bar = sim.spme_outputs.cep_bar(:);
    else
        cep_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cs_ave_n')
        cs_ave_n = sim.spme_outputs.cs_ave_n(:);
    else
        cs_ave_n = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cs_ave_p')
        cs_ave_p = sim.spme_outputs.cs_ave_p(:);
    else
        cs_ave_p = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'c_en_dot_bar')
        c_en_dot_bar = sim.spme_outputs.c_en_dot_bar(:);
    else
        c_en_dot_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'c_es_dot_bar')
        c_es_dot_bar = sim.spme_outputs.c_es_dot_bar(:);
    else
        c_es_dot_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'c_ep_dot_bar')
        c_ep_dot_bar = sim.spme_outputs.c_ep_dot_bar(:);
    else
        c_ep_dot_bar = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cs_ave_n_dot')
        cs_ave_n_dot = sim.spme_outputs.cs_ave_n_dot(:);
    else
        cs_ave_n_dot = nan(size(sim.time));
    end
    if isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'cs_ave_p_dot')
        cs_ave_p_dot = sim.spme_outputs.cs_ave_p_dot(:);
    else
        cs_ave_p_dot = nan(size(sim.time));
    end
    if isfield(sim, 'p') && isstruct(sim.p)
        p = sim.p;
    else
        error('Missing sim.p parameters required for volume-weighted inventory calculation.');
    end

    % Volume-weighted electrolyte inventory [mol] per region and total.
    vol_e_n = p.epsilon_e_n * p.Area_n * p.L_n;
    vol_e_s = p.epsilon_e_s * p.Area_s * p.L_s;
    vol_e_p = p.epsilon_e_p * p.Area_p * p.L_p;
    n_e_n = cen_bar * vol_e_n;
    n_e_s = ces_bar * vol_e_s;
    n_e_p = cep_bar * vol_e_p;
    n_e_total = n_e_n + n_e_s + n_e_p;
    vol_s_n = p.epsilon_s_n * p.Area_n * p.L_n;
    vol_s_p = p.epsilon_s_p * p.Area_p * p.L_p;
    dt = p.delta_t;

    dN_e_n_step = c_en_dot_bar * vol_e_n * dt;
    dN_e_s_step = c_es_dot_bar * vol_e_s * dt;
    dN_e_p_step = c_ep_dot_bar * vol_e_p * dt;
    dN_s_n_step = cs_ave_n_dot * vol_s_n * dt;
    dN_s_p_step = cs_ave_p_dot * vol_s_p * dt;
    dN_total_step = dN_e_n_step + dN_e_s_step + dN_e_p_step + dN_s_n_step + dN_s_p_step;

    figure('Name', 'Sequence signals and cycle capacities', 'Color', 'w');
    tl = tiledlayout(4,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(1);
    plot(sim.time, sim.current, 'm-', 'LineWidth', 1.1); grid on;
    ylabel('Current [A]');
    title('{0.1C discharge, 1C charge, UDDS, 0.1C charge} from SoC = 100%');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(3);
    plot(sim.time, sim.voltage, 'k-', 'LineWidth', 1.1); grid on;
    ylabel('Voltage [V]');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(5);
    plot(sim.time, sim.soc, 'b-', 'LineWidth', 1.1); grid on;
    ylabel('SoC [-]');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(7);
    plot(sim.time, delta_sei * 1e9, 'r-', 'LineWidth', 1.1); grid on;
    ylabel('\delta_{sei} [nm]');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(2);
    plot(sim.time, dN_e_n_step, 'b-', 'LineWidth', 0.9); hold on;
    plot(sim.time, dN_e_s_step, 'm-', 'LineWidth', 0.9);
    plot(sim.time, dN_e_p_step, 'r-', 'LineWidth', 0.9);
    plot(sim.time, dN_s_n_step, 'c-', 'LineWidth', 0.9);
    plot(sim.time, dN_s_p_step, 'y-', 'LineWidth', 0.9);
    plot(sim.time, dN_total_step, 'k-', 'LineWidth', 1.2);
    grid on;
    ylabel('\DeltaN_{Li} per step [mol]');
    legend('\DeltaN_{e,n}', '\DeltaN_{e,s}', '\DeltaN_{e,p}', ...
        '\DeltaN_{s,n}', '\DeltaN_{s,p}', '\DeltaN_{total}', 'Location', 'best');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(4);
    plot(sim.time, cs_ave_n, 'b-', 'LineWidth', 1.0); hold on;
    plot(sim.time, cs_ave_p, 'r-', 'LineWidth', 1.0);
    grid on;
    ylabel('c_{s,ave} [mol/m^3]');
    legend('cs\_ave\_n', 'cs\_ave\_p', 'Location', 'best');
    mark_capacity_cycle_edges(cycle_switch_t);

    nexttile(6);
    plot(sim.time, n_e_n, 'b-', 'LineWidth', 1.0); hold on;
    plot(sim.time, n_e_s, 'm-', 'LineWidth', 1.0);
    plot(sim.time, n_e_p, 'r-', 'LineWidth', 1.0);
    plot(sim.time, n_e_total, 'k-', 'LineWidth', 1.3);
    grid on;
    ylabel('N_e [mol]');
    legend('N_{e,n}', 'N_{e,s}', 'N_{e,p}', 'N_{e,total}', 'Location', 'best');
    mark_capacity_cycle_edges(cycle_switch_t);

    [evt_x, evt_mode, evt_pair_id, evt_label, evt_ah, evt_q_anode, evt_q_cathode, ...
        evt_dNs_n, evt_dNs_p, evt_dNe, evt_li_residual_mol, evt_li_residual_Ah] = ...
        build_capacity_events(sim, cs_ave_n, cs_ave_p, cen_bar, ces_bar, cep_bar, ...
            vol_s_n, vol_s_p, vol_e_n, vol_e_s, vol_e_p, p.Faraday);
    evt_charge_ah = nan(size(evt_x));
    evt_discharge_ah = nan(size(evt_x));
    evt_charge_q_anode = nan(size(evt_x));
    evt_discharge_q_anode = nan(size(evt_x));
    evt_charge_q_cathode = nan(size(evt_x));
    evt_discharge_q_cathode = nan(size(evt_x));
    is_charge = strcmp(evt_mode, 'charge');
    is_discharge = strcmp(evt_mode, 'discharge');
    evt_charge_ah(is_charge) = evt_ah(is_charge);
    evt_discharge_ah(is_discharge) = evt_ah(is_discharge);
    evt_charge_q_anode(is_charge) = evt_q_anode(is_charge);
    evt_discharge_q_anode(is_discharge) = evt_q_anode(is_discharge);
    evt_charge_q_cathode(is_charge) = evt_q_cathode(is_charge);
    evt_discharge_q_cathode(is_discharge) = evt_q_cathode(is_discharge);

    nexttile(8);
    plot(evt_x, evt_charge_ah, 'o-b', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
    plot(evt_x, evt_discharge_ah, 's-r', 'LineWidth', 1.2, 'MarkerSize', 6);
    plot(evt_x, evt_charge_q_anode, '^--c', 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(evt_x, evt_discharge_q_anode, 'v--c', 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(evt_x, evt_charge_q_cathode, '^--m', 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(evt_x, evt_discharge_q_cathode, 'v--m', 'LineWidth', 1.0, 'MarkerSize', 5);
    grid on;
    xticks(evt_x);
    xticklabels(evt_label);
    xtickangle(45);
    xlabel('Capacity event order in profile');
    ylabel('Capacity [Ah]');
    legend('Charged Ah (integral)', 'Discharged Ah (integral)', ...
        'Charged Ah (c_s window, anode)', 'Discharged Ah (c_s window, anode)', ...
        'Charged Ah (c_s window, cathode)', 'Discharged Ah (c_s window, cathode)', ...
        'Location', 'best');
    title(tl, '{0.1C discharge, 1C charge, UDDS, 0.1C charge} from SoC = 100%');

    out = struct();
    out.sim = sim;
    out.cycle_id = cycle_id;
    out.charged_Ah = charged_ah;
    out.discharged_Ah = discharged_ah;
    out.net_Ah = net_ah;
    out.expected_cycles = expected_cycles;
    out.capacity_event_x = evt_x(:);
    out.capacity_event_mode = evt_mode(:);
    out.capacity_event_pair_id = evt_pair_id(:);
    out.capacity_event_label = evt_label(:);
    out.capacity_event_Ah = evt_ah(:);
    out.capacity_event_Q_anode_Ah = evt_q_anode(:);
    out.capacity_event_Q_cathode_Ah = evt_q_cathode(:);
    out.capacity_event_dNs_n_mol = evt_dNs_n(:);
    out.capacity_event_dNs_p_mol = evt_dNs_p(:);
    out.capacity_event_dNe_mol = evt_dNe(:);
    out.capacity_event_li_residual_mol = evt_li_residual_mol(:);
    out.capacity_event_li_residual_Ah = evt_li_residual_Ah(:);
    out.delta_sei = delta_sei(:);
    out.cen_bar = cen_bar(:);
    out.ces_bar = ces_bar(:);
    out.cep_bar = cep_bar(:);
    out.cs_ave_n = cs_ave_n(:);
    out.cs_ave_p = cs_ave_p(:);
    out.n_e_n = n_e_n(:);
    out.n_e_s = n_e_s(:);
    out.n_e_p = n_e_p(:);
    out.n_e_total = n_e_total(:);
    out.vol_e_n = vol_e_n;
    out.vol_e_s = vol_e_s;
    out.vol_e_p = vol_e_p;
    out.vol_s_n = vol_s_n;
    out.vol_s_p = vol_s_p;
    out.dN_e_n_step = dN_e_n_step(:);
    out.dN_e_s_step = dN_e_s_step(:);
    out.dN_e_p_step = dN_e_p_step(:);
    out.dN_s_n_step = dN_s_n_step(:);
    out.dN_s_p_step = dN_s_p_step(:);
    out.dN_total_step = dN_total_step(:);
    out.degradation_model = degradation_model;
end

function mark_capacity_cycle_edges(cycle_switch_t)
    for i = 1:numel(cycle_switch_t)
        xline(cycle_switch_t(i), '--', sprintf('cycle %d', i+1), ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
end

function [evt_x, evt_mode, evt_pair_id, evt_label, evt_ah, evt_q_anode, evt_q_cathode, ...
        evt_dNs_n, evt_dNs_p, evt_dNe, evt_li_residual_mol, evt_li_residual_Ah] = ...
        build_capacity_events(sim, cs_ave_n, cs_ave_p, cen_bar, ces_bar, cep_bar, ...
            vol_s_n, vol_s_p, vol_e_n, vol_e_s, vol_e_p, faraday_const)
    seg_log = sim.segment_log;
    seg_idx_all = sim.segment_index(:);
    n_seg = numel(seg_log);

    evt_x = [];
    evt_mode = {};
    evt_pair_id = [];
    evt_label = {};
    evt_ah = [];
    evt_q_anode = [];
    evt_q_cathode = [];
    evt_dNs_n = [];
    evt_dNs_p = [];
    evt_dNe = [];
    evt_li_residual_mol = [];
    evt_li_residual_Ah = [];

    evt_counter = 0;
    first_mode = '';
    pair_counter = 1;
    for s = 1:n_seg
        if strcmp(seg_log(s).name, 'REST')
            continue;
        end
        idx_start = find(seg_idx_all == s, 1, 'first');
        idx_end_profile = find(seg_idx_all == s, 1, 'last');
        if isempty(idx_start) || isempty(idx_end_profile)
            continue;
        end

        idx_end_window = idx_end_profile;
        if s < n_seg && strcmp(seg_log(s+1).name, 'REST')
            idx_end_rest = find(seg_idx_all == (s+1), 1, 'last');
            if ~isempty(idx_end_rest)
                idx_end_window = idx_end_rest;
            end
        end

        if strcmp(seg_log(s).name, 'CC')
            is_charge = seg_log(s).crate < 0;
        else
            % Dynamic profiles are treated as discharge-oriented.
            is_charge = false;
        end

        dNs_n_signed = (cs_ave_n(idx_end_window) - cs_ave_n(idx_start)) * vol_s_n;
        dNs_p_signed = (cs_ave_p(idx_end_window) - cs_ave_p(idx_start)) * vol_s_p;
        dNe_signed = ((cen_bar(idx_end_window) - cen_bar(idx_start)) * vol_e_n) + ...
            ((ces_bar(idx_end_window) - ces_bar(idx_start)) * vol_e_s) + ...
            ((cep_bar(idx_end_window) - cep_bar(idx_start)) * vol_e_p);
        li_residual_mol = dNs_n_signed + dNs_p_signed + dNe_signed;
        li_residual_Ah = li_residual_mol * faraday_const / 3600;

        if is_charge
            mode = 'charge';
            ah_val = seg_log(s).charged_Ah;
            dn_n = cs_ave_n(idx_end_window) - cs_ave_n(idx_start);
            dn_p = cs_ave_p(idx_start) - cs_ave_p(idx_end_window);
        else
            mode = 'discharge';
            ah_val = seg_log(s).discharged_Ah;
            dn_n = cs_ave_n(idx_start) - cs_ave_n(idx_end_window);
            dn_p = cs_ave_p(idx_end_window) - cs_ave_p(idx_start);
        end

        evt_counter = evt_counter + 1;
        if isempty(first_mode)
            first_mode = mode;
            pair_id = pair_counter;
        elseif strcmp(mode, first_mode)
            pair_counter = pair_counter + 1;
            pair_id = pair_counter;
        else
            pair_id = pair_counter;
        end

        if strcmp(mode, 'charge')
            mode_tag = 'cha';
        else
            mode_tag = 'discha';
        end
        evt_x(evt_counter, 1) = evt_counter; %#ok<AGROW>
        evt_mode{evt_counter, 1} = mode; %#ok<AGROW>
        evt_pair_id(evt_counter, 1) = pair_id; %#ok<AGROW>
        evt_label{evt_counter, 1} = sprintf('%d,%s', pair_id, mode_tag); %#ok<AGROW>
        evt_ah(evt_counter, 1) = ah_val; %#ok<AGROW>
        evt_q_anode(evt_counter, 1) = max(dn_n, 0) * vol_s_n * faraday_const / 3600; %#ok<AGROW>
        evt_q_cathode(evt_counter, 1) = max(dn_p, 0) * vol_s_p * faraday_const / 3600; %#ok<AGROW>
        evt_dNs_n(evt_counter, 1) = dNs_n_signed; %#ok<AGROW>
        evt_dNs_p(evt_counter, 1) = dNs_p_signed; %#ok<AGROW>
        evt_dNe(evt_counter, 1) = dNe_signed; %#ok<AGROW>
        evt_li_residual_mol(evt_counter, 1) = li_residual_mol; %#ok<AGROW>
        evt_li_residual_Ah(evt_counter, 1) = li_residual_Ah; %#ok<AGROW>
    end
end

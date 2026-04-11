function out = test_sim_profile_sequence_fuds_cc3()
%TEST_SIM_PROFILE_SEQUENCE_FUDS_CC3
% Validate sequence-repeat feature with 3 cycles of:
%   FUDS profile -> 1C CC charge
% and visualize current, voltage, SoC, and SEI thickness.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    % Conditions
    model_name = 'SPMe';
    init_soc = 0.4;
    sequence = {'FUDS', 'CC', -1};
    n_cycles = 3;
    t_amb = 25;
    delta_t = 1;
    Nr = 20; Nxn = 20; Nxs = 10; Nxp = 20;
    deg_model = 'ramadass';

    sim = sim_profile_sequence_LGM50(model_name, init_soc, sequence, ...
        'Repeat', n_cycles, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', deg_model, 'SaveOutput', false);

    if strcmpi(deg_model, 'none')
        delta_sei = nan(size(sim.time));
    else
        delta_sei = sim.x(end, :).';
    end

    cycle_edges = find(diff(sim.segment_cycle) > 0) + 1;
    cycle_switch_t = sim.time(cycle_edges);

    figure('Name', 'Sequence repeat test: (FUDS -> 1C charge) x3', 'Color', 'w');

    subplot(4,1,1);
    plot(sim.time, sim.current, 'm-', 'LineWidth', 1.1); grid on;
    ylabel('Current [A]');
    title('(FUDS -> 1C charge) repeated 3 cycles');
    mark_cycle_edges(cycle_switch_t);

    subplot(4,1,2);
    plot(sim.time, sim.voltage, 'k-', 'LineWidth', 1.1); grid on;
    ylabel('Voltage [V]');
    mark_cycle_edges(cycle_switch_t);

    subplot(4,1,3);
    plot(sim.time, sim.soc, 'b-', 'LineWidth', 1.1); grid on;
    ylabel('SoC [-]');
    mark_cycle_edges(cycle_switch_t);

    subplot(4,1,4);
    plot(sim.time, delta_sei * 1e9, 'r-', 'LineWidth', 1.1); grid on;
    xlabel('Time [s]');
    ylabel('\delta_{sei} [nm]');
    mark_cycle_edges(cycle_switch_t);

    out = struct();
    out.sim = sim;
    out.time = sim.time(:);
    out.current = sim.current(:);
    out.voltage = sim.voltage(:);
    out.soc = sim.soc(:);
    out.delta_sei = delta_sei(:);
    out.cycle_switch_time = cycle_switch_t(:);
    out.conditions = struct('model', model_name, ...
                            'init_soc', init_soc, ...
                            'sequence', {sequence}, ...
                            'repeat', n_cycles, ...
                            'degradation_model', deg_model, ...
                            'Tamb_degC', t_amb, ...
                            'DeltaT_s', delta_t, ...
                            'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp);
end

function mark_cycle_edges(cycle_switch_t)
    for i = 1:numel(cycle_switch_t)
        xline(cycle_switch_t(i), '--', sprintf('cycle %d', i+1), ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
end

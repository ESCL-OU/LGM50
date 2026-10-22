function out = test_spme_sei()
%TEST_SPME_SEI
% SPMe SEI trajectory aligned with legacy ode_spmet SEI test:
% 1C charge (SOC 0 -> Vmax), then 1C discharge (-> Vmin), one continuous run.

    tests_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(tests_dir);
    addpath(repo_dir);
    addpath(fullfile(repo_dir, 'SPMe'));
    addpath(fullfile(repo_dir, 'param'));

    disp('Single Particle Model w/ Electrolyte (SPMe)')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % Match legacy-comparison settings.
    t_amb = 25;
    delta_t = 1;
    init_soc = 0.0;
    CRate_charge = -1.0;
    CRate_discharge = 1.0;
    Nr = 20; Nxn = 20; Nxs = 10; Nxp = 20;

    sim = sim_profile_sequence_LGM50('SPMe', init_soc, {'CC', CRate_charge, 'CC', CRate_discharge}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'ramadass', ...
        'SPMeOutputs', {'j_sei', 'delta_sei'}, ...
        'SaveOutput', false);

    nT = numel(sim.time);
    if ~isfield(sim, 'spme_outputs')
        error('sim.spme_outputs missing; request SPMeOutputs in simulation call.');
    end
    jsei = sim.spme_outputs.j_sei(:);
    dsei = sim.spme_outputs.delta_sei(:);

    k_switch = find(sim.segment_index(:) == 2, 1, 'first');
    if isempty(k_switch)
        k_switch = find(diff(sign(sim.current(:))) ~= 0, 1, 'first') + 1;
    end
    if isempty(k_switch) || k_switch < 1 || k_switch > nT
        k_switch = nT;
    end

    figure('Name', 'LGM50-main SPMe SEI (charge -> discharge)', 'Color', 'w');
    subplot(4,1,1);
    plot(sim.time, sim.voltage, 'k-', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    title('SPMe trajectory (charge -> discharge)');
    xline(sim.time(k_switch), '--', 'Phase switch');

    subplot(4,1,2);
    plot(sim.time, sim.current, 'm-', 'LineWidth', 1.2); grid on;
    ylabel('Current [A]');
    xline(sim.time(k_switch), '--');
    title('Applied current');

    subplot(4,1,3);
    plot(sim.time, jsei, 'b-', 'LineWidth', 1.2); grid on;
    ylabel('j_{sei} [mol/(m^2 s)]');
    xline(sim.time(k_switch), '--');
    title('SEI side-reaction flux');

    subplot(4,1,4);
    plot(sim.time, dsei*1e9, 'r-', 'LineWidth', 1.2); grid on;
    xlabel('Time [s]');
    ylabel('\delta_{sei} [nm]');
    xline(sim.time(k_switch), '--');
    title('SEI thickness');

    out = struct();
    out.time = sim.time(:);
    out.voltage = sim.voltage(:);
    out.current = sim.current(:);
    out.j_sei = jsei(:);
    out.delta_sei = dsei(:);
    out.sim = sim;
    out.info = struct('repo', 'LGM50-main', 'model', 'SPMe', ...
                      'degradation_model', 'ramadass', ...
                      'init_soc', init_soc, ...
                      'CRate_charge', CRate_charge, ...
                      'CRate_discharge', CRate_discharge, ...
                      'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
                      'delta_t', delta_t, ...
                      'k_switch', k_switch, ...
                      't_switch', sim.time(k_switch));
end

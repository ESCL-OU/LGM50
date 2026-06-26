function out = test_dfn_degradation_coupling_comparison()
% Compare DFN degradation behavior between spme_compatible and yang_total_flux.

degradation_model = 'safari';
enable_lithium_plating = true;
enable_anode_porosity_loss = true;

init_soc = 0.1;
sequence1 = {'REST', 0.25, 'CC', -0.2, 'REST', 0.25, 'CC', 0.2};
sequence2 = {'CC', -1.0, 'CC', 1.0};
%block_sequence = [sequence1, repmat(sequence2, 1, 20)];
block_sequence = {'CC', -1.0, 'CC', 1.0};
delta_t = 1;
t_amb = 25;
Nr = 10;
Nxn = 10;
Nxs = 5;
Nxp = 10;
modes = {'spme_compatible', 'yang_total_flux'};
out = struct();

for m = 1:numel(modes)
    mode = modes{m};

    sim = sim_profile_sequence_LGM50('DFN', init_soc, block_sequence, ...
        'Repeat', 1, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', degradation_model, ...
        'EnableLithiumPlating', enable_lithium_plating, ...
        'EnableAnodePorosityLoss', enable_anode_porosity_loss, ...
        'DFNLplCouplingMode', mode, ...
        'SaveOutput', false, ...
        'ShowProgress', true, ...
        'MemoryLight', false);

    fprintf('\nMode requested: %s\n', mode);
    fprintf('Mode stored in sim.p.dfn_lpl_coupling_mode: %s\n', sim.p.dfn_lpl_coupling_mode);
    if isfield(sim, 'dfn_lpl_coupling_mode')
        fprintf('Mode stored in sim.dfn_lpl_coupling_mode: %s\n', sim.dfn_lpl_coupling_mode);
    end

    out.(mode).sim = sim;
end

plot_dfn_mode_comparison(out, modes);

end

function plot_dfn_mode_comparison(out, modes)

figure('Name', 'DFN coupling mode comparison', 'Color', 'w');
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for m = 1:numel(modes)
    mode = modes{m};
    sim = out.(mode).sim;
    t_hr = sim.time(:) / 3600;
    d = sim.dfn_outputs;

    disp(struct('mode',mode,'time',numel(t_hr),'LDP_sep',numel(d.LDP_sep_n),'j_lpl_sep',numel(d.j_lpl_sep_n)))

    nexttile(1); hold on;
    plot(t_hr, sim.voltage(:), 'DisplayName', mode);
    ylabel('Voltage [V]'); grid on;

    nexttile(2); hold on;
    plot(t_hr, sim.soc(:), 'DisplayName', mode);
    ylabel('SOC [-]'); grid on;

    nexttile(3); hold on;
    y_LDP = d.LDP_sep_n(:);
    n = min(numel(t_hr), numel(y_LDP));
    plot(t_hr(1:n), y_LDP(1:n), 'DisplayName', mode);
    ylabel('Separator-side LDP [V]'); grid on;

    nexttile(4); hold on;
    y_j = d.j_lpl_sep_n(:);
    n = min(numel(t_hr), numel(y_j));
    plot(t_hr(1:n), y_j(1:n), 'DisplayName', mode);    
    ylabel('Separator-side j_{lpl} [mol/m^2/s]'); grid on;

    nexttile(5); hold on;
    y_dli = d.delta_li_sep_n(:);
    n = min(numel(t_hr), numel(y_dli));
    plot(t_hr(1:n), y_dli(1:n) * 1e9, 'DisplayName', mode);
    ylabel('\delta_{Li,sep} [nm]'); xlabel('Time [h]'); grid on;

    nexttile(6); hold on;
    y_eps = d.epsilon_e_sep_n(:);
    n = min(numel(t_hr), numel(y_eps));
    plot(t_hr(1:n), y_eps(1:n), 'DisplayName', mode);
    ylabel('\epsilon_{e,n,sep} [-]'); xlabel('Time [h]'); grid on;
end

for k = 1:6
    nexttile(k);
    legend('Location', 'best', 'Interpreter', 'none');
end

figure('Name', 'DFN final spatial profiles', 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for m = 1:numel(modes)
    mode = modes{m};
    sim = out.(mode).sim;
    d = sim.dfn_outputs;
    x_n = linspace(0, 1, size(d.LDP_n, 1));

    nexttile(1); hold on;
    plot(x_n, d.LDP_n(:, end), 'DisplayName', mode);
    ylabel('LDP_n [V]'); grid on;

    nexttile(2); hold on;
    plot(x_n, d.j_lpl_n(:, end), 'DisplayName', mode);
    ylabel('j_{lpl,n} [mol/m^2/s]'); grid on;

    nexttile(3); hold on;
    plot(x_n, d.delta_li_n(:, end) * 1e9, 'DisplayName', mode);
    ylabel('\delta_{Li,n} [nm]'); xlabel('x/L_n'); grid on;

    nexttile(4); hold on;
    plot(x_n, d.epsilon_e_n_dyn(:, end), 'DisplayName', mode);
    ylabel('\epsilon_{e,n} [-]'); xlabel('x/L_n'); grid on;
end

for k = 1:4
    nexttile(k);
    legend('Location', 'best', 'Interpreter', 'none');
end

end
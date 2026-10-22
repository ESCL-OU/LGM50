function out = test_electrolyte_mass_balance()
%TEST_ELECTROLYTE_MASS_BALANCE
% Investigate electrolyte lithium mass balance during 1C charge from SoC 0%.
%
% This test runs SPMe with CC charge, then audits the electrolyte mass
% balance decomposition at each time step using audit_electrolyte_mass_balance.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    model_name = 'SPMe';
    init_soc = 0.0;
    c_rate = -1; % sign convention in this repo: negative = charge
    delta_t = 1;
    t_amb = 25;

    sim = sim_profile_LGM50(model_name, init_soc, {'CC', c_rate}, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'DegradationModel', 'none', ...
        'SPMeOutputs', {'jn', 'jp'}, ...
        'SaveOutput', false);

    assert(isfield(sim, 'c_e') && ~isempty(sim.c_e), ...
        'sim.c_e is required for electrolyte mass-balance audit.');
    assert(isfield(sim, 'spme_outputs') && isfield(sim.spme_outputs, 'jn') ...
        && isfield(sim.spme_outputs, 'jp'), ...
        'sim.spme_outputs.jn/jp are required for electrolyte mass-balance audit.');

    t_vec = sim.time(:);
    c_e_hist = sim.c_e;
    jn_hist = sim.spme_outputs.jn(:);
    jp_hist = sim.spme_outputs.jp(:);
    p = sim.p;

    Nt = numel(t_vec);
    dNe_total = zeros(Nt,1);
    dNe_t1 = zeros(Nt,1);
    dNe_t2 = zeros(Nt,1);
    dNe_src = zeros(Nt,1);
    dNe_n = zeros(Nt,1);
    dNe_s = zeros(Nt,1);
    dNe_p = zeros(Nt,1);

    for k = 1:Nt
        A = audit_electrolyte_mass_balance(c_e_hist(:,k), jn_hist(k), jp_hist(k), p);
        dNe_total(k) = A.dNe_total;
        dNe_t1(k) = A.terms.term1_dD_grad2;
        dNe_t2(k) = A.terms.term2_D_lap;
        dNe_src(k) = A.terms.source;
        dNe_n(k) = A.regions.anode.total;
        dNe_s(k) = A.regions.sep.total;
        dNe_p(k) = A.regions.cathode.total;
    end

    figure('Name', 'Electrolyte mass-balance audit: 1C charge from SoC 0%', 'Color', 'w');
    tl = tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t_vec, dNe_total, 'k', 'LineWidth', 1.4); hold on;
    plot(t_vec, dNe_t1, '--', 'LineWidth', 1.1);
    plot(t_vec, dNe_t2, '--', 'LineWidth', 1.1);
    plot(t_vec, dNe_src, '--', 'LineWidth', 1.1);
    grid on;
    ylabel('dN_e/dt [mol/m^2/s]');
    legend('total', 'term1: dD(dc/dx)^2', 'term2: D d2c/dx2', 'source', ...
        'Location', 'best');
    title('Global electrolyte inventory rate decomposition');

    nexttile;
    plot(t_vec, dNe_n, 'b-', 'LineWidth', 1.1); hold on;
    plot(t_vec, dNe_s, 'm-', 'LineWidth', 1.1);
    plot(t_vec, dNe_p, 'r-', 'LineWidth', 1.1);
    plot(t_vec, dNe_n + dNe_s + dNe_p, 'k--', 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Regional dN_e/dt [mol/m^2/s]');
    legend('anode', 'separator', 'cathode', 'sum(regions)', 'Location', 'best');
    title('Regional electrolyte inventory rates');
    title(tl, 'test\_electrolyte\_mass\_balance');

    out = struct();
    out.sim = sim;
    out.time = t_vec;
    out.dNe_total = dNe_total;
    out.dNe_term1 = dNe_t1;
    out.dNe_term2 = dNe_t2;
    out.dNe_source = dNe_src;
    out.dNe_anode = dNe_n;
    out.dNe_separator = dNe_s;
    out.dNe_cathode = dNe_p;
    out.conditions = struct('model', model_name, 'init_soc', init_soc, ...
        'c_rate', c_rate, 'DeltaT_s', delta_t, 'Tamb_degC', t_amb);
end

function out = test_spme_yang_vs_safari_charge(assert_alignment)
%TEST_SPME_YANG_VS_SAFARI_CHARGE Compare Yang and Safari SEI variables during charging.
%
% Usage:
%   out = test_spme_yang_vs_safari_charge();
%   out = test_spme_yang_vs_safari_charge(false);  % compare only, no asserts

    if nargin < 1 || isempty(assert_alignment)
        assert_alignment = true;
    end

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    % Charge-only protocol
    init_soc = 0.1;
    c_rate_charge = -0.5;
    t_amb = 25;
    delta_t = 1;
    Nr = 15; Nxn = 30; Nxs = 5; Nxp = 30;

    % Keep Safari settings explicit for reproducibility.
    nsei = 10;
    epsilon_sei = 1.0;

    req_out = {'j_sei', 'j_side', 'delta_sei', 'delta_sei_dot', ...
               'c_s_ec', 'i0s', 'eta_n_side', 'phi_se', 'eta_s'};

    fprintf('[%s] Running Yang charge simulation...\n', datestr(now, 31));
    sim_yang = sim_profile_LGM50('SPMe', init_soc, {'CC', c_rate_charge}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'yang', ...
        'ElectrolyteModel', 'fdm', ...
        'SPMeOutputs', req_out, ...
        'SaveOutput', false, ...
        'SuppressSummary', true);

    fprintf('[%s] Running Safari charge simulation...\n', datestr(now, 31));
    sim_safari = sim_profile_LGM50('SPMe', init_soc, {'CC', c_rate_charge}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'safari', ...
        'ElectrolyteModel', 'fdm', ...
        'Nsei', nsei, 'EpsilonSEI', epsilon_sei, ...
        'SPMeOutputs', req_out, ...
        'SaveOutput', false, ...
        'SuppressSummary', true);

    vars_to_compare = {'j_sei', 'j_side', 'delta_sei', 'delta_sei_dot', ...
                       'c_s_ec', 'i0s', 'eta_n_side', 'phi_se', 'eta_s'};
    ylabels = {'j_{sei} [mol/(m^2 s)]', 'j_{side} [mol/(m^2 s)]', ...
               '\delta_{sei} [m]', 'd\delta_{sei}/dt [m/s]', ...
               'c_{s,ec} [mol/m^3]', 'i_{0s} [A/m^2]', ...
               '\eta_n^{side} [V]', '\phi_{se} [V]', '\eta_s [V]'};

    % Relative RMSE thresholds against Yang baseline.
    tol.rel_rmse.j_sei = 0.25;
    tol.rel_rmse.j_side = 0.25;
    tol.rel_rmse.delta_sei = 0.25;
    tol.rel_rmse.delta_sei_dot = 0.30;
    tol.rel_rmse.c_s_ec = 0.15;
    tol.rel_rmse.i0s = 0.15;
    tol.rel_rmse.eta_n_side = 0.05;
    tol.rel_rmse.phi_se = 0.05;
    tol.rel_rmse.eta_s = 0.05;

    metrics = struct();
    nvars = numel(vars_to_compare);
    figure('Name', 'SEI variables: Yang vs Safari (charge only)', 'Color', 'w');
    tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:nvars
        vn = vars_to_compare{i};
        y_yang = sim_yang.spme_outputs.(vn)(:);
        y_safari = sim_safari.spme_outputs.(vn)(:);

        [t_cmp, y_ref, y_cmp] = overlap_series(sim_yang.time(:), y_yang, sim_safari.time(:), y_safari);
        rmse_abs = sqrt(mean((y_ref - y_cmp).^2, 'omitnan'));
        ref_scale = max(abs(y_ref));
        ref_scale = max(ref_scale, eps);
        rmse_rel = rmse_abs / ref_scale;

        end_ref = y_ref(end);
        end_cmp = y_cmp(end);
        if abs(end_ref) < eps
            end_ratio = NaN;
        else
            end_ratio = end_cmp / end_ref;
        end

        metrics.(vn) = struct();
        metrics.(vn).rmse_abs = rmse_abs;
        metrics.(vn).rmse_rel = rmse_rel;
        metrics.(vn).end_yang = end_ref;
        metrics.(vn).end_safari = end_cmp;
        metrics.(vn).end_ratio_safari_over_yang = end_ratio;
        metrics.(vn).time = t_cmp;
        metrics.(vn).yang = y_ref;
        metrics.(vn).safari = y_cmp;

        nexttile;
        plot(t_cmp, y_ref, 'b-', 'LineWidth', 1.3); hold on;
        plot(t_cmp, y_cmp, 'r--', 'LineWidth', 1.2); grid on;
        title(sprintf('%s (rel RMSE=%.3g)', vn, rmse_rel), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel(ylabels{i});
        if i == 1
            legend('Yang', 'Safari', 'Location', 'best');
        end
    end

    fprintf('\nYang vs Safari (charge-only) SEI comparison metrics:\n');
    fprintf('----------------------------------------------------\n');
    for i = 1:nvars
        vn = vars_to_compare{i};
        m = metrics.(vn);
        fprintf('%-14s | relRMSE=%9.3e | absRMSE=%9.3e | end_ratio(Safari/Yang)=%9.3e\n', ...
            vn, m.rmse_rel, m.rmse_abs, m.end_ratio_safari_over_yang);
    end
    fprintf('----------------------------------------------------\n\n');

    if assert_alignment
        bad = {};
        for i = 1:nvars
            vn = vars_to_compare{i};
            limit = tol.rel_rmse.(vn);
            if metrics.(vn).rmse_rel > limit
                bad{end+1} = sprintf('%s relRMSE %.3e > %.3e', ... %#ok<AGROW>
                    vn, metrics.(vn).rmse_rel, limit);
            end
        end
        if ~isempty(bad)
            error(['Yang vs Safari mismatch exceeded thresholds during charge:\n  %s\n' ...
                   'Set assert_alignment=false to inspect plots/metrics without failing.'], ...
                   strjoin(bad, '\n  '));
        end
    end

    out = struct();
    out.sim_yang = sim_yang;
    out.sim_safari = sim_safari;
    out.metrics = metrics;
    out.thresholds = tol;
    out.settings = struct('init_soc', init_soc, ...
                          'c_rate_charge', c_rate_charge, ...
                          't_amb', t_amb, ...
                          'delta_t', delta_t, ...
                          'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
                          'nsei', nsei, 'epsilon_sei', epsilon_sei);
end

function [t_cmp, y1_cmp, y2_cmp] = overlap_series(t1, y1, t2, y2)
    t_end = min(t1(end), t2(end));
    t_cmp = t1(t1 <= t_end);
    y1_cmp = y1(1:numel(t_cmp));
    y2_cmp = interp1(t2, y2, t_cmp, 'linear');
end

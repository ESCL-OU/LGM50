function out = debug_spme_dfca_values()
%DEBUG_SPME_DFCA_VALUES
% Evaluate dfca_n/dfca_s/dfca_p along a CC simulation trajectory.
% Also estimates how much the dfca terms contribute to V_electrolytePolar.

    import casadi.*

    base_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'param'));

    % Match baseline settings in legacy-vs-current test.
    init_SOC = 0.0;
    CRate = -1.0;
    t_amb = 25;
    delta_t = 1;
    Nr = 30; Nxn = 70; Nxs = 35; Nxp = 70;
    disable_thermal = true;

    p = load_SPMe_params(init_SOC, t_amb, Nr, delta_t, Nxn, Nxs, Nxp, disable_thermal, 'none');
    I = p.OneC * CRate;
    NT = max(2, floor(1.5 * 3600 / abs(CRate)) + 1);
    t = (0:NT-1)' * delta_t;

    m = SPMe(p, 'idas');
    x = p.x0;

    dfca_n = zeros(NT,1);
    dfca_s = zeros(NT,1);
    dfca_p = zeros(NT,1);
    vpol_no_dfca = zeros(NT,1);
    vpol_with_dfca = zeros(NT,1);
    dvpol_dfca_only = zeros(NT,1);
    v_total = zeros(NT,1);

    t_last = NT;
    for k = 1:NT
        % Use current state x to evaluate c_ex and voltage components.
        cex_struct = m.f_cex('x', x);
        c_ex = full(cex_struct.c_ex);
        y_struct = m.h('x', x, 'p', I);
        y = full(y_struct.y);
        v_total(k) = y(1);

        [dfca_n(k), dfca_s(k), dfca_p(k), vpol_no_dfca(k), vpol_with_dfca(k)] = local_dfca_and_vpol(c_ex, p);
        dvpol_dfca_only(k) = vpol_with_dfca(k) - vpol_no_dfca(k);

        if k < NT
            [x, y_next] = m.step(x, I);
            if y_next(1) >= p.volt_max
                t_last = k;
                break;
            end
        end
    end

    % Trim to simulated horizon.
    t = t(1:t_last);
    dfca_n = dfca_n(1:t_last);
    dfca_s = dfca_s(1:t_last);
    dfca_p = dfca_p(1:t_last);
    vpol_no_dfca = vpol_no_dfca(1:t_last);
    vpol_with_dfca = vpol_with_dfca(1:t_last);
    dvpol_dfca_only = dvpol_dfca_only(1:t_last);
    v_total = v_total(1:t_last);

    % Numeric summary table for quick inspection.
    summary = table( ...
        [dfca_n(1); mean(dfca_n); max(dfca_n)], ...
        [dfca_s(1); mean(dfca_s); max(dfca_s)], ...
        [dfca_p(1); mean(dfca_p); max(dfca_p)], ...
        [dvpol_dfca_only(1); mean(dvpol_dfca_only); max(dvpol_dfca_only)], ...
        'VariableNames', {'dfca_n','dfca_s','dfca_p','DeltaVpol_dfcaOnly'}, ...
        'RowNames', {'initial','mean','max'});
    disp(summary);

    figure('Name', 'dfca diagnostics in SPMe', 'Color', 'w');
    subplot(4,1,1);
    plot(t, dfca_n, 'b-', 'LineWidth', 1.2); hold on;
    plot(t, dfca_s, 'k--', 'LineWidth', 1.1);
    plot(t, dfca_p, 'r-.', 'LineWidth', 1.1); grid on;
    ylabel('dfca [-]');
    title('Electrolyte activity coefficients');
    legend('dfca_n','dfca_s','dfca_p','Location','best');

    subplot(4,1,2);
    plot(t, vpol_no_dfca, 'k--', 'LineWidth', 1.1); hold on;
    plot(t, vpol_with_dfca, 'm-', 'LineWidth', 1.2); grid on;
    ylabel('V_{elytePolar} [V]');
    title('Polarization term with and without dfca');
    legend('without dfca','with dfca','Location','best');

    subplot(4,1,3);
    plot(t, dvpol_dfca_only, 'g-', 'LineWidth', 1.2); grid on;
    ylabel('\DeltaV_{polar} [V]');
    title('Contribution from dfca terms only');

    subplot(4,1,4);
    plot(t, v_total, 'b-', 'LineWidth', 1.2); grid on;
    xlabel('Time [s]');
    ylabel('V [V]');
    title('Total cell voltage');

    out = struct();
    out.time = t;
    out.dfca_n = dfca_n;
    out.dfca_s = dfca_s;
    out.dfca_p = dfca_p;
    out.vpol_no_dfca = vpol_no_dfca;
    out.vpol_with_dfca = vpol_with_dfca;
    out.dvpol_dfca_only = dvpol_dfca_only;
    out.v_total = v_total;
    out.summary = summary;
end

function [dfn, dfs, dfp, vpol0, vpol1] = local_dfca_and_vpol(c_ex, p)
    ce0n = c_ex(1);
    cens = c_ex(p.Nxn+1);
    cesp = c_ex(p.Nxn+p.Nxs+1);
    ce0p = c_ex(end);

    cen_bar = mean(c_ex(1:p.Nxn+1));
    ces_bar = mean(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1)));
    cep_bar = mean(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1)));

    dfn = electrolyteAct(cen_bar, p.T_amb, p);
    dfs = electrolyteAct(ces_bar, p.T_amb, p);
    dfp = electrolyteAct(cep_bar, p.T_amb, p);

    base = (2*p.R*p.T_amb)/(p.Faraday);
    dlog_n = log(cens) - log(ce0n);
    dlog_s = log(cesp) - log(cens);
    dlog_p = log(ce0p) - log(cesp);

    vpol0 = base * (dlog_n + dlog_s + dlog_p);
    vpol1 = base * ((1+dfn)*dlog_n + (1+dfs)*dlog_s + (1+dfp)*dlog_p);
end

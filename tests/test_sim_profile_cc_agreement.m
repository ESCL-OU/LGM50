function out = test_sim_profile_cc_agreement()
%TEST_SIM_PROFILE_CC_AGREEMENT
% Compare sim_profile_LGM50('CC') against legacy CC implementations:
%   - sim_SPMe_CC.m style loop
%   - sim_DFN_CC.m style loop
%   - sim_ECM_CC.m style loop

    import casadi.*

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'DFN'));
    addpath(fullfile(base_dir, 'ECM'));
    addpath(fullfile(base_dir, 'param'));

    CRate = -1;
    t_amb = 25;
    init_SOC = 0.0;
    delta_t = 1;
    disable_thermal = true;

    sim_spme_profile = sim_profile_LGM50('SPMe', init_SOC, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, 'SaveOutput', false);
    sim_dfn_profile = sim_profile_LGM50('DFN', init_SOC, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, 'SaveOutput', false);
    sim_ecm_profile = sim_profile_LGM50('ECM', init_SOC, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, 'SaveOutput', false);

    [t_spme_ref, v_spme_ref, soc_spme_ref] = run_spme_cc_reference(init_SOC, CRate, t_amb, delta_t, disable_thermal);
    [t_dfn_ref,  v_dfn_ref,  soc_dfn_ref ] = run_dfn_cc_reference(init_SOC, CRate, t_amb, delta_t, disable_thermal);
    [t_ecm_ref,  v_ecm_ref,  soc_ecm_ref ] = run_ecm_cc_reference(init_SOC, CRate, t_amb, delta_t);

    [rmse_spme, t_spme_cmp, v_spme_profile_cmp, v_spme_ref_cmp] = overlap_rmse(sim_spme_profile.time, sim_spme_profile.voltage, t_spme_ref, v_spme_ref);
    [rmse_dfn,  t_dfn_cmp,  v_dfn_profile_cmp,  v_dfn_ref_cmp ] = overlap_rmse(sim_dfn_profile.time,  sim_dfn_profile.voltage,  t_dfn_ref,  v_dfn_ref);
    [rmse_ecm,  t_ecm_cmp,  v_ecm_profile_cmp,  v_ecm_ref_cmp ] = overlap_rmse(sim_ecm_profile.time,  sim_ecm_profile.voltage,  t_ecm_ref,  v_ecm_ref);
    [~, ~, soc_spme_profile_cmp, soc_spme_ref_cmp] = overlap_rmse(sim_spme_profile.time, sim_spme_profile.soc, t_spme_ref, soc_spme_ref);
    [~, ~, soc_dfn_profile_cmp,  soc_dfn_ref_cmp ] = overlap_rmse(sim_dfn_profile.time,  sim_dfn_profile.soc,  t_dfn_ref,  soc_dfn_ref);
    [~, ~, soc_ecm_profile_cmp,  soc_ecm_ref_cmp ] = overlap_rmse(sim_ecm_profile.time,  sim_ecm_profile.soc,  t_ecm_ref,  soc_ecm_ref);

    figure('Name', 'sim\_profile\_LGM50 vs legacy CC', 'Color', 'w');

    subplot(3,2,1);
    plot(t_spme_cmp, v_spme_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_spme_cmp, v_spme_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    title(sprintf('SPMe CC charge (CRate=%g), RMSE = %.5f V', CRate, rmse_spme));
    legend('sim\_profile\_LGM50', 'sim\_SPMe\_CC-style', 'Location', 'best');

    subplot(3,2,2);
    plot(t_spme_cmp, soc_spme_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_spme_cmp, soc_spme_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('SoC [-]');
    title('SPMe SoC');
    legend('sim\_profile\_LGM50', 'sim\_SPMe\_CC-style', 'Location', 'best');

    subplot(3,2,3);
    plot(t_dfn_cmp, v_dfn_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_dfn_cmp, v_dfn_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    title(sprintf('DFN CC charge (CRate=%g), RMSE = %.5f V', CRate, rmse_dfn));
    legend('sim\_profile\_LGM50', 'sim\_DFN\_CC-style', 'Location', 'best');

    subplot(3,2,4);
    plot(t_dfn_cmp, soc_dfn_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_dfn_cmp, soc_dfn_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('SoC [-]');
    title('DFN SoC');
    legend('sim\_profile\_LGM50', 'sim\_DFN\_CC-style', 'Location', 'best');

    subplot(3,2,5);
    plot(t_ecm_cmp, v_ecm_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_ecm_cmp, v_ecm_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    xlabel('Time [s]');
    title(sprintf('ECM CC charge (CRate=%g), RMSE = %.5f V', CRate, rmse_ecm));
    legend('sim\_profile\_LGM50', 'sim\_ECM\_CC-style', 'Location', 'best');

    subplot(3,2,6);
    plot(t_ecm_cmp, soc_ecm_profile_cmp, 'b-', 'LineWidth', 1.4); hold on;
    plot(t_ecm_cmp, soc_ecm_ref_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('SoC [-]');
    xlabel('Time [s]');
    title('ECM SoC');
    legend('sim\_profile\_LGM50', 'sim\_ECM\_CC-style', 'Location', 'best');

    out = struct();
    out.conditions = struct('CRate', CRate, 't_amb', t_amb, 'init_SOC', init_SOC, 'delta_t', delta_t);
    out.rmse = struct('SPMe', rmse_spme, 'DFN', rmse_dfn, 'ECM', rmse_ecm);
    out.profile = struct('SPMe', sim_spme_profile, 'DFN', sim_dfn_profile, 'ECM', sim_ecm_profile);
    out.reference = struct( ...
        'SPMe', struct('time', t_spme_ref, 'voltage', v_spme_ref, 'soc', soc_spme_ref), ...
        'DFN',  struct('time', t_dfn_ref,  'voltage', v_dfn_ref,  'soc', soc_dfn_ref), ...
        'ECM',  struct('time', t_ecm_ref,  'voltage', v_ecm_ref,  'soc', soc_ecm_ref));
end

function [t_ref, v_ref, soc_ref] = run_spme_cc_reference(init_SOC, CRate, t_amb, delta_t, disable_thermal)
    sim_ref = sim_profile_LGM50('SPMe', init_SOC, {'CC', CRate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'DegradationModel', 'none', 'SaveOutput', false);
    t_ref = sim_ref.time(:);
    v_ref = sim_ref.voltage(:);
    soc_ref = sim_ref.soc(:);
end

function [t_ref, v_ref, soc_ref] = run_dfn_cc_reference(init_SOC, CRate, t_amb, delta_t, disable_thermal)
    p = load_DFN_params(init_SOC, t_amb, delta_t, disable_thermal);
    t_ref = (0:delta_t:floor(3600/abs(CRate))).';
    NT = numel(t_ref);
    I = p.OneC * CRate * ones(NT,1);

    x = zeros(length(p.x0), NT);
    z = zeros(length(p.z0), NT);
    y = zeros(length(p.y0), NT);
    x(:,1) = p.x0; z(:,1) = p.z0; y(:,1) = p.y0;

    t_last = NT;
    for k = 1:NT-1
        if k == 1
            Cur_vec = [I(k), I(k), I(k+1)];
        else
            Cur_vec = [I(k-1), I(k), I(k+1)];
        end
        [x(:,k+1), z(:,k+1), y(:,k+1)] = step_dfn(x(:,k), z(:,k), Cur_vec, p, 0);

        Vk = y(1, k+1);
        if Vk <= p.volt_min || Vk >= p.volt_max
            t_last = k+1;
            break;
        end
    end
    t_ref = t_ref(1:t_last);
    v_ref = y(1,1:t_last).';
    soc_ref = y(2,1:t_last).';
end

function [t_ref, v_ref, soc_ref] = run_ecm_cc_reference(init_SOC, CRate, t_amb, delta_t)
    n_order = 2;
    [p, ecm_p, ocv_p] = load_ECM_params(init_SOC, n_order);
    p.delta_t = delta_t;
    p.T_amb = t_amb + 273.15;
    ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');

    t_ref = (0:delta_t:floor(3600/abs(CRate))).';
    NT = numel(t_ref);
    I = p.OneC * CRate;

    x = zeros(length(p.x0), NT);
    v_ref = zeros(NT, 1);
    x(:,1) = p.x0;
    v_ref(1) = p.y0;

    for k = 1:NT-1
        [x(:,k+1), v_ref(k+1)] = ecm_ss.step(x(:,k), I);
    end
    soc_ref = x(1,:).';
end

function [rmse_v, t_cmp, v1_cmp, v2_cmp] = overlap_rmse(t1, v1, t2, v2)
    t_end = min(t1(end), t2(end));
    t_cmp = t1(t1 <= t_end);
    v1_cmp = v1(1:numel(t_cmp));
    v2_cmp = interp1(t2, v2, t_cmp, 'linear');
    rmse_v = sqrt(mean((v1_cmp - v2_cmp).^2, 'omitnan'));
end

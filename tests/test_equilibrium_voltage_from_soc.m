function result = test_equilibrium_voltage_from_soc(target_soc)
%TEST_EQUILIBRIUM_VOLTAGE_FROM_SOC Validate single-SoC equilibrium voltage.

    if nargin < 1
        target_soc = 0.5;
    end
    if ~isscalar(target_soc)
        error('target_soc must be a single scalar value.');
    end
    if target_soc < 0 || target_soc > 1
        error('target_soc must be in [0, 1].');
    end

    calc = calculate_equilibrium_voltage_from_soc(target_soc);
    p = calc.p;

    soc_from_anode = (calc.cn_interp - calc.cn_soc0) / (calc.cn_soc100 - calc.cn_soc0);
    soc_from_cathode = (calc.cp_interp - calc.cp_soc0) / (calc.cp_soc100 - calc.cp_soc0);
    soc_interp_diff_pct = 100 * (soc_from_anode - soc_from_cathode);

    p_check = p;
    p_check.csn0 = calc.cn_interp;
    p_check.csp0 = calc.cp_interp;
    [cn_back, cp_back] = init_cs(p_check, calc.Veq);

    cn_err = cn_back - calc.cn_interp;
    cp_err = cp_back - calc.cp_interp;
    cn_err_pct = 100 * cn_err / max(abs(calc.cn_interp), eps);
    cp_err_pct = 100 * cp_err / max(abs(calc.cp_interp), eps);

    result = struct();
    result.SOC_target = calc.SOC_target;
    result.Veq = calc.Veq;
    result.cn_interp = calc.cn_interp;
    result.cp_interp = calc.cp_interp;
    result.soc_from_anode_interp = soc_from_anode;
    result.soc_from_cathode_interp = soc_from_cathode;
    result.soc_interp_diff_pct = soc_interp_diff_pct;
    result.cn_initcs = cn_back;
    result.cp_initcs = cp_back;
    result.cn_err = cn_err;
    result.cp_err = cp_err;
    result.cn_err_pct = cn_err_pct;
    result.cp_err_pct = cp_err_pct;

    fprintf('\n=== Single-SoC Equilibrium Voltage Test ===\n');
    fprintf('Target SoC: %.4f\n', calc.SOC_target);
    fprintf('Veq from interpolation: %.6f V\n', calc.Veq);
    fprintf('SOC from anode interpolation  : %.6f\n', soc_from_anode);
    fprintf('SOC from cathode interpolation: %.6f\n', soc_from_cathode);
    fprintf('SOC interpolation diff        : %.3e %%\n', soc_interp_diff_pct);
    fprintf('cn interp / init_cs / err     : %.6f / %.6f / %.3e (%.3e %%)\n', ...
        calc.cn_interp, cn_back, cn_err, cn_err_pct);
    fprintf('cp interp / init_cs / err     : %.6f / %.6f / %.3e (%.3e %%)\n\n', ...
        calc.cp_interp, cp_back, cp_err, cp_err_pct);
end

function results = test_init_cs_init_models()
%TEST_INIT_CS_INIT_MODELS Legacy-name test for SoC-based initialization.

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    param_dir = fullfile(base_dir, 'param');
    addpath(base_dir);

    if ~isfolder(param_dir)
        error('param folder not found relative to %s', base_dir);
    end

    addpath(param_dir);

    soc_list = [0.1, 0.15, 0.2, 0.3, 0.4];
    n = numel(soc_list);
    soc_col = zeros(n, 1);
    v_eq_col = zeros(n, 1);
    csn_ref_col = zeros(n, 1);
    csp_ref_col = zeros(n, 1);
    csn_est_col = zeros(n, 1);
    csp_est_col = zeros(n, 1);

    for k = 1:n
        soc = soc_list(k);
        calc = calculate_equilibrium_voltage_from_soc(soc);
        p = calc.p;
        [csn_est, csp_est] = init_cs(p, calc.Veq);

        soc_col(k) = soc;
        v_eq_col(k) = calc.Veq;
        csn_ref_col(k) = calc.cn_interp;
        csp_ref_col(k) = calc.cp_interp;
        csn_est_col(k) = csn_est;
        csp_est_col(k) = csp_est;
    end

    csn_abs_err = csn_est_col - csn_ref_col;
    csp_abs_err = csp_est_col - csp_ref_col;
    csn_rel_err = 100 * csn_abs_err ./ csn_ref_col;
    csp_rel_err = 100 * csp_abs_err ./ csp_ref_col;

    results = table( ...
        soc_col, v_eq_col, ...
        csn_ref_col, csn_est_col, csn_abs_err, csn_rel_err, ...
        csp_ref_col, csp_est_col, csp_abs_err, csp_rel_err, ...
        'VariableNames', { ...
            'SOC', 'Veq', ...
            'csn_ref', 'csn_init_cs', 'csn_abs_err', 'csn_rel_err_pct', ...
            'csp_ref', 'csp_init_cs', 'csp_abs_err', 'csp_rel_err_pct'});

    disp('init_cs consistency check vs SoC-based calculate_equilibrium_voltage_from_soc:');
    disp(results);
end

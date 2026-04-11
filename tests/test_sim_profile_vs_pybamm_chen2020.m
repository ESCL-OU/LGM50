function out = test_sim_profile_vs_pybamm_chen2020()
%TEST_SIM_PROFILE_VS_PYBAMM_CHEN2020
% Compare MATLAB implementation (sim_profile_LGM50) vs PyBaMM (Chen2020)
% for 1C CC operation from a chosen SoC, using SPMe and DFN.
%
% Requirements:
%   - Python environment configured in MATLAB
%   - pybamm package installed in that Python environment

    tests_dir = fileparts(mfilename('fullpath'));
    base_dir = fileparts(tests_dir);
    addpath(base_dir);
    addpath(fullfile(base_dir, 'SPMe'));
    addpath(fullfile(base_dir, 'DFN'));
    addpath(fullfile(base_dir, 'param'));
    ensure_python_runtime_unloaded();

    init_soc = 1.0;
    c_rate = 1; % sign convention in this repo: positive=discharge, negative=charge
    t_amb = 25;
    delta_t = 1;

    python_exe = resolve_python_executable();
    [pybamm_ok, py_diag] = check_pybamm_available_external(python_exe);
    if ~pybamm_ok
        error(['PyBaMM is not available for external Python execution.\n', ...
               'Python executable: %s\n', ...
               'Probe output:\n%s\n', ...
               'Install with:\n', ...
               '  "%s" -m pip install pybamm'], ...
              py_diag.python_executable, py_diag.probe_output, py_diag.python_executable);
    end

    % MATLAB implementation results
    sim_spme = sim_profile_LGM50('SPMe', init_soc, {'CC', c_rate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'DegradationModel', 'none', ...
        'ElectrolyteModel', 'fvm', ...
        'SaveOutput', false);
    sim_dfn = sim_profile_LGM50('DFN', init_soc, {'CC', c_rate}, ...
        'Tamb', t_amb, 'DeltaT', delta_t, ...
        'SaveOutput', false);

    % PyBaMM (Chen2020) results
    py_spme = run_pybamm_cc_chen2020_external('SPMe', init_soc, c_rate, python_exe);
    py_dfn = run_pybamm_cc_chen2020_external('DFN', init_soc, c_rate, python_exe);

    [rmse_v_spme, t_spme_cmp, v_spme_mat_cmp, v_spme_py_cmp] = ...
        overlap_rmse(sim_spme.time, sim_spme.voltage, py_spme.time_s, py_spme.voltage_v);
    [rmse_v_dfn, t_dfn_cmp, v_dfn_mat_cmp, v_dfn_py_cmp] = ...
        overlap_rmse(sim_dfn.time, sim_dfn.voltage, py_dfn.time_s, py_dfn.voltage_v);
    [rmse_soc_spme, ~, soc_spme_mat_cmp, soc_spme_py_cmp] = ...
        overlap_rmse(sim_spme.time, sim_spme.soc, py_spme.time_s, py_spme.soc);
    [rmse_soc_dfn, ~, soc_dfn_mat_cmp, soc_dfn_py_cmp] = ...
        overlap_rmse(sim_dfn.time, sim_dfn.soc, py_dfn.time_s, py_dfn.soc);

    if c_rate >= 0
        mode_text = 'discharge';
    else
        mode_text = 'charge';
    end
    figure('Name', sprintf('MATLAB vs PyBaMM (Chen2020), %.3gC %s from SoC=%.2f', ...
        abs(c_rate), mode_text, init_soc), 'Color', 'w');

    subplot(2,2,1);
    plot(t_spme_cmp, v_spme_mat_cmp, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_spme_cmp, v_spme_py_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    title(sprintf('SPMe Voltage (RMSE = %.5f V)', rmse_v_spme));
    legend('sim\_profile\_LGM50', 'PyBaMM Chen2020', 'Location', 'best');

    subplot(2,2,2);
    plot(t_spme_cmp, soc_spme_mat_cmp, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_spme_cmp, soc_spme_py_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('SoC [-]');
    title(sprintf('SPMe SoC (RMSE = %.5f)', rmse_soc_spme));
    legend('sim\_profile\_LGM50', 'PyBaMM Chen2020', 'Location', 'best');

    subplot(2,2,3);
    plot(t_dfn_cmp, v_dfn_mat_cmp, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_dfn_cmp, v_dfn_py_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('Voltage [V]');
    xlabel('Time [s]');
    title(sprintf('DFN Voltage (RMSE = %.5f V)', rmse_v_dfn));
    legend('sim\_profile\_LGM50', 'PyBaMM Chen2020', 'Location', 'best');

    subplot(2,2,4);
    plot(t_dfn_cmp, soc_dfn_mat_cmp, 'b-', 'LineWidth', 1.3); hold on;
    plot(t_dfn_cmp, soc_dfn_py_cmp, 'r--', 'LineWidth', 1.2); grid on;
    ylabel('SoC [-]');
    xlabel('Time [s]');
    title(sprintf('DFN SoC (RMSE = %.5f)', rmse_soc_dfn));
    legend('sim\_profile\_LGM50', 'PyBaMM Chen2020', 'Location', 'best');

    out = struct();
    out.conditions = struct('init_soc', init_soc, 'c_rate', c_rate, ...
        'mode', mode_text, ...
        'Tamb_degC', t_amb, 'DeltaT_s', delta_t);
    out.rmse = struct( ...
        'SPMe_voltage', rmse_v_spme, ...
        'SPMe_soc', rmse_soc_spme, ...
        'DFN_voltage', rmse_v_dfn, ...
        'DFN_soc', rmse_soc_dfn);
    out.matlab = struct('SPMe', sim_spme, 'DFN', sim_dfn);
    out.pybamm = struct('SPMe', py_spme, 'DFN', py_dfn);
end

function py_out = run_pybamm_cc_chen2020_external(model_name, init_soc, c_rate, python_exe)
    script_path = [tempname, '.py'];
    out_path = [tempname, '.json'];
    cleanup_obj = onCleanup(@() cleanup_temp_files({script_path, out_path})); %#ok<NASGU>

    write_pybamm_runner_script(script_path);
    cmd = sprintf('"%s" "%s" "%s" "%.16g" "%.16g" "%s"', ...
        python_exe, script_path, upper(char(model_name)), init_soc, c_rate, out_path);
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('PyBaMM external run failed for %s:\n%s', model_name, cmdout);
    end
    if ~isfile(out_path)
        error('PyBaMM external run did not produce output file: %s', out_path);
    end

    data = jsondecode(fileread(out_path));
    py_out = struct();
    py_out.time_s = data.time_s(:);
    py_out.voltage_v = data.voltage_v(:);
    py_out.soc = data.soc(:);
end

function [rmse_v, t_cmp, v1_cmp, v2_cmp] = overlap_rmse(t1, v1, t2, v2)
    t_end = min(t1(end), t2(end));
    t_cmp = t1(t1 <= t_end);
    v1_cmp = v1(1:numel(t_cmp));
    v2_cmp = interp1(t2, v2, t_cmp, 'linear');
    rmse_v = sqrt(mean((v1_cmp - v2_cmp).^2, 'omitnan'));
end

function [ok, diag] = check_pybamm_available_external(python_exe)
    ok = false;
    diag = struct('python_executable', python_exe, 'probe_output', '');
    cmd = sprintf('"%s" -c "import pybamm; print(pybamm.__version__)"', python_exe);
    [status, cmdout] = system(cmd);
    diag.probe_output = strtrim(cmdout);
    ok = (status == 0);
end

function python_exe = resolve_python_executable()
    preferred = '/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/bin/python3';
    python_exe = preferred;
    if isfile(preferred)
        return;
    end
    try
        pe = pyenv;
        if isprop(pe, 'Executable') && ~isempty(pe.Executable) && isfile(char(pe.Executable))
            python_exe = char(pe.Executable);
            return;
        end
    catch
        % keep fallback
    end
    python_exe = 'python3';
end

function write_pybamm_runner_script(script_path)
    fid = fopen(script_path, 'w');
    if fid < 0
        error('Could not create temporary Python runner script: %s', script_path);
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = {
        'import os'
        'import sys'
        'import json'
        'os.environ.setdefault("PYBAMM_DISABLE_TELEMETRY", "true")'
        'os.environ.setdefault("PYBAMM_ENABLE_TELEMETRY", "0")'
        'import pybamm'
        ''
        'model_name = sys.argv[1].upper()'
        'init_soc = float(sys.argv[2])'
        'c_rate = float(sys.argv[3])'
        'c_rate_abs = abs(c_rate)'
        'out_path = sys.argv[4]'
        ''
        'if model_name == "SPME":'
        '    model = pybamm.lithium_ion.SPMe()'
        'elif model_name == "DFN":'
        '    model = pybamm.lithium_ion.DFN()'
        'else:'
        '    raise ValueError(f"Unsupported model_name: {model_name}")'
        ''
        'param = pybamm.ParameterValues("Chen2020")'
        'try:'
        '    param = param.set_initial_stoichiometries(init_soc)'
        'except Exception:'
        '    pass'
        ''
        'if c_rate < 0:'
        '    step = f"Charge at {c_rate_abs} C until 4.2 V"'
        'elif c_rate > 0:'
        '    step = f"Discharge at {c_rate_abs} C until 2.5 V"'
        'else:'
        '    raise ValueError("c_rate must be non-zero")'
        'experiment = pybamm.Experiment([step])'
        'sim = pybamm.Simulation(model, parameter_values=param, experiment=experiment)'
        'try:'
        '    sol = sim.solve(initial_soc=init_soc)'
        'except Exception:'
        '    sol = sim.solve()'
        ''
        'if not hasattr(sol, "__getitem__"):'
        '    raise RuntimeError("PyBaMM returned EmptySolution (likely infeasible step for init_soc/protocol).")'
        't = sol["Time [s]"].entries'
        'v = sol["Terminal voltage [V]"].entries'
        ''
        'soc = None'
        'for name in ["SoC", "SOC", "State of Charge", "State of charge"]:'
        '    try:'
        '        soc = sol[name].entries'
        '        break'
        '    except Exception:'
        '        pass'
        ''
        'if soc is None:'
        '    q_dis = sol["Discharge capacity [A.h]"].entries'
        '    try:'
        '        q_nom = float(param["Nominal cell capacity [A.h]"])'
        '    except Exception:'
        '        q_nom = 5.0'
        '    soc = init_soc - q_dis / q_nom'
        ''
        'soc = [float(min(1.0, max(0.0, s))) for s in soc]'
        'data = {'
        '    "time_s": [float(x) for x in t],'
        '    "voltage_v": [float(x) for x in v],'
        '    "soc": soc'
        '}'
        'with open(out_path, "w", encoding="utf-8") as f:'
        '    json.dump(data, f)'
    };
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
end

function cleanup_temp_files(paths)
    for i = 1:numel(paths)
        p = paths{i};
        if ischar(p) || isstring(p)
            if isfile(p)
                delete(p);
            end
        end
    end
end

function ensure_python_runtime_unloaded()
    % Avoid CasADi shared-library collision between MATLAB CasADi MEX and
    % Python casadi loaded via MATLAB's in-process Python runtime.
    try
        pe = pyenv;
        if isprop(pe, 'Status') && strcmpi(char(pe.Status), 'Loaded')
            try
                terminate(pe);
            catch ME
                error(['MATLAB Python runtime is loaded and could not be terminated.\n', ...
                       'This can conflict with MATLAB CasADi MEX.\n', ...
                       'Please restart MATLAB, then run this test before any py.* calls.\n', ...
                       'Original terminate error: %s'], ME.message);
            end
        end
    catch
        % If pyenv is unavailable, continue.
    end
end

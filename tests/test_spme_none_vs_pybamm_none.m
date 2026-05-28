%% test_spme_none_vs_pybamm_none.m
% Compare your SPMe.m (no SEI) vs PyBaMM SPMe (SEI = none).
%
% Requirements:
%   - MATLAB model files on path
%   - Python environment with pybamm installed

clear; clc; close all;

%% User settings
C_rate = -1;              % + discharge, - charge in this repo
init_soc = 0.10;
dt = 1;

Nr = 20;
Nxn = 20;
Nxs = 10;
Nxp = 20;
Tamb_C = 25;

pybamm_param_set = "Chen2020";

%% Paths
tests_dir = fileparts(mfilename("fullpath"));
project_root = fileparts(tests_dir);
addpath(tests_dir);
addpath(project_root);

subdirs = {"param", "SPMe", "DFN", "ECM", "ECM 2"};
for k = 1:numel(subdirs)
    candidate = fullfile(project_root, subdirs{k});
    if isfolder(candidate)
        addpath(candidate);
    end
end

%% Build initial SPMe parameters (no SEI model)
disable_thermal = true;
degradation_model = "none";
electrolyte_model = "fdm";

p = load_SPMe_params(init_soc, Tamb_C, Nr, dt, Nxn, Nxs, Nxp, ...
    disable_thermal, degradation_model, electrolyte_model);

%% Run your SPMe model (no SEI)
sim_my = sim_profile_LGM50("SPMe", init_soc, {"CC", C_rate}, ...
    "Nxn", Nxn, ...
    "Nxs", Nxs, ...
    "Nxp", Nxp, ...
    "Nr", Nr, ...
    "DeltaT", dt, ...
    "Tamb", Tamb_C, ...
    "DegradationModel", "none", ...
    "ElectrolyteModel", "fdm", ...
    "UseElectrolyteActivity", false, ...
    "SuppressSummary", true);

t_my = sim_my.time;
V_my = sim_my.voltage;
SOC_my = sim_my.soc;
I_my = sim_my.current;
Voc_my = zeros(size(t_my));
Vrxn_my = zeros(size(t_my));
Vohm_e_my = zeros(size(t_my));
Vconc_my = zeros(size(t_my));
V_recon_my = zeros(size(t_my));

spme_pot = SPMe(sim_my.p, "idas");
for k = 1:numel(t_my)
    [eta_n_k, eta_p_k, Unref_k, Upref_k, ~, Vohm_e_k, Vconc_k] = ...
        spme_pot.potentials(sim_my.x(:, k), I_my(k));
    Voc_my(k) = Upref_k - Unref_k;
    Vrxn_my(k) = eta_p_k - eta_n_k;
    Vohm_e_my(k) = Vohm_e_k;
    Vconc_my(k) = Vconc_k;
    V_recon_my(k) = Voc_my(k) + Vrxn_my(k) + Vohm_e_my(k) + Vconc_my(k);
end

%% Run PyBaMM SPMe with SEI disabled
pybamm = py.importlib.import_module("pybamm");
np = py.importlib.import_module("numpy");

options = py.dict(pyargs( ...
    "SEI", "none", ...
    "thermal", "isothermal"));

model = pybamm.lithium_ion.SPMe(options);
param = pybamm.ParameterValues(pybamm_param_set);
try_update_param(param, {"Nominal cell capacity [A.h]"}, p.OneC);
try_update_param(param, {"Typical current [A]"}, p.OneC);

%% Parameter comparison (capacity-driving terms)
F = 96485.33212;
Q_spme_from_soc_window = F * p.epsilon_s_n * p.L_n * p.Area_n * ...
    (p.c_n_soc100 - p.c_n_soc0) / 3600;
Q_pybamm_nom = get_pybamm_param(param, "Nominal cell capacity [A.h]", NaN);
Vmin_pybamm = get_pybamm_param(param, "Lower voltage cut-off [V]", NaN);
Vmax_pybamm = get_pybamm_param(param, "Upper voltage cut-off [V]", NaN);
csn_max_pybamm = get_pybamm_param(param, "Maximum concentration in negative electrode [mol.m-3]", NaN);
csp_max_pybamm = get_pybamm_param(param, "Maximum concentration in positive electrode [mol.m-3]", NaN);

fprintf("\nCapacity-related parameter comparison:\n");
fprintf("  SPMe OneC                       = %.6g Ah\n", p.OneC);
fprintf("  SPMe capacity from SOC window   = %.6g Ah\n", Q_spme_from_soc_window);
fprintf("  SPMe voltage cutoffs [min, max] = [%.6g, %.6g] V\n", p.volt_min, p.volt_max);
fprintf("  SPMe c_s_max [n, p]             = [%.6g, %.6g] mol/m3\n", p.c_s_n_max, p.c_s_p_max);
fprintf("  PyBaMM nominal capacity         = %.6g Ah\n", Q_pybamm_nom);
fprintf("  PyBaMM voltage cutoffs [min,max]= [%.6g, %.6g] V\n", Vmin_pybamm, Vmax_pybamm);
fprintf("  PyBaMM c_s_max [n, p]           = [%.6g, %.6g] mol/m3\n", csn_max_pybamm, csp_max_pybamm);

if C_rate >= 0
    V_end = p.volt_min;
    step = sprintf("Discharge at %.6gC until %.6g V", abs(C_rate), V_end);
else
    V_end = p.volt_max;
    step = sprintf("Charge at %.6gC until %.6g V", abs(C_rate), V_end);
end

experiment = pybamm.Experiment({step}, pyargs("period", sprintf("%.6g seconds", dt)));
sim_py = pybamm.Simulation(model, ...
    pyargs("parameter_values", param, "experiment", experiment));
sol = sim_py.solve(pyargs("initial_soc", init_soc));

t_py = double(py.array.array("d", py.numpy.nditer(sol.t)));
V_py = get_pybamm_var(sol, { ...
    "Terminal voltage [V]", ...
    "Voltage [V]"});
I_py = get_pybamm_var(sol, { ...
    "Current [A]", ...
    "Current variable [A]", ...
    "Cell current [A]"});
Voc_py = get_pybamm_var(sol, { ...
    "Surface open-circuit voltage [V]", ...
    "Battery open-circuit voltage [V]", ...
    "Bulk open-circuit voltage [V]"});
Vrxn_py = get_pybamm_var(sol, { ...
    "X-averaged reaction overpotential [V]", ...
    "X-averaged battery reaction overpotential [V]"});
Vohm_e_py = get_pybamm_var(sol, { ...
    "X-averaged electrolyte ohmic losses [V]", ...
    "X-averaged battery electrolyte ohmic losses [V]"});
Vconc_py = get_pybamm_var(sol, { ...
    "X-averaged concentration overpotential [V]", ...
    "X-averaged battery concentration overpotential [V]"});
V_recon_py = Voc_py + Vrxn_py + Vohm_e_py + Vconc_py;

SOC_py = get_pybamm_var(sol, { ...
    "SoC", ...
    "SOC", ...
    "State of Charge", ...
    "State of charge"});
if isempty(SOC_py)
    fprintf("PyBaMM SOC variable not found; using fallback from Discharge capacity [A.h].\n");
    q_dis = get_pybamm_var(sol, {"Discharge capacity [A.h]"});
    try
        q_nom = double(param.get("Nominal cell capacity [A.h]"));
        fprintf("Using PyBaMM parameter for fallback SOC: Nominal cell capacity [A.h] = %.6g\n", q_nom);
    catch
        q_nom = p.OneC;
        fprintf("Fallback Nominal capacity not found in PyBaMM; using p.OneC = %.6g Ah\n", q_nom);
    end
    SOC_py = init_soc - q_dis / q_nom;
end

%% Numerical comparison
V_py_i = interp1(t_py, V_py, t_my, "linear", "extrap");
SOC_py_i = interp1(t_py, SOC_py, t_my, "linear", "extrap");
I_py_i = interp1(t_py, I_py, t_my, "linear", "extrap");

fprintf("\nFinal values:\n");
fprintf("  My SPMe voltage     = %.6f V\n", V_my(end));
fprintf("  PyBaMM voltage      = %.6f V\n", V_py(end));
fprintf("  My SPMe SOC         = %.6f\n", SOC_my(end));
fprintf("  PyBaMM SOC          = %.6f\n", SOC_py(end));
fprintf("  My SPMe current     = %.6f A\n", I_my(end));
fprintf("  PyBaMM current      = %.6f A\n", I_py(end));

fprintf("\nRMSE:\n");
fprintf("  Voltage = %.6e V\n", sqrt(mean((V_my - V_py_i).^2)));
fprintf("  SOC     = %.6e\n", sqrt(mean((SOC_my - SOC_py_i).^2)));
fprintf("  Current = %.6e A\n", sqrt(mean((I_my - I_py_i).^2)));

fprintf("\nTime horizon:\n");
fprintf("  MATLAB SPMe end time = %.2f s\n", t_my(end));
fprintf("  PyBaMM end time      = %.2f s\n", t_py(end));
if isfield(sim_my, "terminated_early") && sim_my.terminated_early
    fprintf("  MATLAB termination reason: %s\n", sim_my.termination_reason);
end

%% Plots
figure;
plot(t_my/3600, V_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, V_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("Voltage [V]");
title(sprintf("SPMe none-SEI voltage comparison, %.2gC", C_rate));
legend("My SPMe none-SEI", "PyBaMM none-SEI", "Location", "best");

figure;
plot(t_my/3600, I_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, I_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("Current [A]");
title(sprintf("SPMe none-SEI current comparison, %.2gC", C_rate));
legend("My SPMe none-SEI", "PyBaMM none-SEI", "Location", "best");

figure;
plot(t_my/3600, SOC_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, SOC_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SOC [-]");
title(sprintf("SPMe none-SEI SOC comparison, %.2gC", C_rate));
legend("My SPMe none-SEI", "PyBaMM none-SEI", "Location", "best");

%% Voltage decomposition comparison
figure("Name", sprintf("Voltage decomposition, %.2gC", C_rate), "Color", "w");

subplot(3,2,1);
plot(t_my/3600, Voc_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, Voc_py, "--", "LineWidth", 1.3);
grid on; ylabel("OCP [V]");
title("Open-circuit term");
legend("SPMe", "PyBaMM", "Location", "best");

subplot(3,2,2);
plot(t_my/3600, Vrxn_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, Vrxn_py, "--", "LineWidth", 1.3);
grid on; ylabel("V_{rxn} [V]");
title("Reaction overpotential");
legend("SPMe", "PyBaMM", "Location", "best");

subplot(3,2,3);
plot(t_my/3600, Vohm_e_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, Vohm_e_py, "--", "LineWidth", 1.3);
grid on; ylabel("V_{ohm,e} [V]");
title("Electrolyte ohmic term");
legend("SPMe", "PyBaMM", "Location", "best");

subplot(3,2,4);
plot(t_my/3600, Vconc_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, Vconc_py, "--", "LineWidth", 1.3);
grid on; ylabel("V_{conc} [V]");
title("Concentration overpotential");
legend("SPMe", "PyBaMM", "Location", "best");

subplot(3,2,5);
plot(t_my/3600, V_recon_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, V_recon_py, "--", "LineWidth", 1.3);
grid on; xlabel("Time [h]"); ylabel("V_{recon} [V]");
title("Reconstructed terminal voltage");
legend("SPMe", "PyBaMM", "Location", "best");

subplot(3,2,6);
plot(t_my/3600, V_my - V_recon_my, "LineWidth", 1.3); hold on;
plot(t_py/3600, V_py - V_recon_py, "--", "LineWidth", 1.3);
grid on; xlabel("Time [h]"); ylabel("Residual [V]");
title("V - V_{recon}");
legend("SPMe", "PyBaMM", "Location", "best");

%% Local function
function try_update_param(param, candidate_keys, val)
    if isstring(candidate_keys)
        key_list = cellstr(candidate_keys);
    elseif iscell(candidate_keys)
        key_list = cell(size(candidate_keys));
        for ii = 1:numel(candidate_keys)
            key_list{ii} = char(string(candidate_keys{ii}));
        end
    else
        key_list = {char(string(candidate_keys))};
    end

    for i = 1:numel(key_list)
        key = key_list{i};
        try
            param.update(py.dict(pyargs(key, val)), pyargs("check_already_exists", true));
            fprintf("Updated PyBaMM parameter: %s = %.6g\n", key, val);
            return;
        catch
        end
    end
    fprintf("Skipped PyBaMM parameter aliases: %s\n", strjoin(key_list, " | "));
end

function y = get_pybamm_var(sol, candidate_names)
    op = py.importlib.import_module("operator");
    np = py.importlib.import_module("numpy");

    for i = 1:numel(candidate_names)
        name = candidate_names{i};
        try
            v = op.getitem(sol, char(name));   % Python: sol[name]
            entries = np.asarray(v.entries).flatten();
            y = double(py.array.array("d", py.numpy.nditer(entries)));
            fprintf("Using PyBaMM variable: %s\n", name);
            return;
        catch ME
            fprintf("Failed PyBaMM variable alias: %s\n", name);
            fprintf("  %s\n", ME.message);
        end
    end

    fprintf("No PyBaMM variable found for aliases: %s\n", strjoin(cellstr(string(candidate_names)), " | "));
    y = [];
end

function val = get_pybamm_param(param, key, fallback)
    try
        val = double(param.get(key));
        fprintf("Using PyBaMM parameter: %s = %.6g\n", key, val);
    catch ME
        fprintf("Failed PyBaMM parameter lookup: %s\n", key);
        fprintf("  %s\n", ME.message);
        val = fallback;
    end
end

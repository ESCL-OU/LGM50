%% test_spme_ramadass_vs_pybamm_reaction_limited.m
% Compare your SPMe.m Ramadass SEI model vs PyBaMM reaction-limited SEI.
%
% Mapping:
%   Your model:       DegradationModel = "ramadass"
%   PyBaMM model:     options("SEI") = "reaction limited"
%
% This follows the same structure as test_spme_yang_vs_pybamm_yang.m.
%
% Requirements:
%   - Your MATLAB model files on path
%   - Python environment with pybamm installed
%
% In MATLAB, if needed:
%   pyenv("Version","/path/to/python")

clear; clc; close all;

%% User settings
C_rate = -1.0;              % + discharge, - charge in your code
init_soc = 0.10;
T_end_s = 3600;
dt = 1;

initial_sei_m = 1.0e-9;     % force both models to same initial SEI

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

%% Build initial SPMe state with matched initial SEI
disable_thermal = true;
degradation_model = "ramadass";
electrolyte_model = "fdm";

p = load_SPMe_params(init_soc, Tamb_C, Nr, dt, Nxn, Nxs, Nxp, ...
    disable_thermal, degradation_model, electrolyte_model);

calc = calculate_equilibrium_voltage_from_soc(init_soc);
[csn0, csp0] = init_cs(p, calc.Veq);

c_n0 = csn0 * ones(p.Nr - 1, 1);
c_p0 = csp0 * ones(p.Nr - 1, 1);
c_e0 = p.c_e * ones(p.Nx - 3, 1);

init_state = struct();
init_state.x = [c_n0; c_p0; c_e0; initial_sei_m];

%% Run your SPMe.m Ramadass model
sim_my = sim_profile_LGM50("SPMe", init_soc, {"CC", C_rate}, ...
    "Nxn", Nxn, ...
    "Nxs", Nxs, ...
    "Nxp", Nxp, ...
    "Nr", Nr, ...
    "DeltaT", dt, ...
    "Tamb", Tamb_C, ...
    "DegradationModel", "ramadass", ...
    "ElectrolyteModel", "fdm", ...
    "UseElectrolyteActivity", false, ...
    "InitialState", init_state, ...
    "SPMeOutputs", {"delta_sei","delta_sei_dot","j_sei","j_side", ...
                    "eta_n_side","phi_se","eta_s","V","SOC"}, ...
    "SuppressSummary", true);

keep = sim_my.time <= T_end_s;

t_my = sim_my.time(keep);
L_my = sim_my.spme_outputs.delta_sei(keep);
dLdt_my = sim_my.spme_outputs.delta_sei_dot(keep);
j_my = sim_my.spme_outputs.j_sei(keep);   % mol/m2/s in your code
V_my = sim_my.voltage(keep);
I_my = sim_my.current(keep);

%% Run PyBaMM from MATLAB Python interface
pybamm = py.importlib.import_module("pybamm");
np = py.importlib.import_module("numpy");

options = py.dict(pyargs( ...
    "SEI", "reaction limited", ...
    "SEI film resistance", "distributed", ...
    "thermal", "isothermal"));

model = pybamm.lithium_ion.SPMe(options);
param = pybamm.ParameterValues(pybamm_param_set);

% Match key SEI parameters to your params_LGM50.m where possible.
% Your Ramadass branch uses:
%   p.i0s      [A/m2]
%   p.Us       [V]
%   p.kappa_P  [S/m] -> PyBaMM resistivity is [Ohm m] = 1/kappa_P
%   p.M_P      [kg/mol]
%   p.rho_P    [kg/m3]
try_update_param(param, {"Initial SEI thickness [m]"}, initial_sei_m);
try_update_param(param, {"SEI resistivity [Ohm.m]"}, 1 / p.kappa_P);
try_update_param(param, {"SEI reaction exchange current density [A.m-2]", ...
                         "SEI reaction exchange current density [A/m2]", ...
                         "SEI exchange current density [A.m-2]"}, p.i0s);
try_update_param(param, {"SEI open-circuit potential [V]", ...
                         "SEI reaction open-circuit potential [V]", ...
                         "SEI equilibrium potential [V]"}, p.Us);
try_update_param(param, {"SEI partial molar volume [m3.mol-1]"}, p.M_P / p.rho_P);

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

%% Extract PyBaMM variables
L_py = get_pybamm_var(sol, { ...
    "X-averaged negative SEI thickness [m]", ...
    "X-averaged negative total SEI thickness [m]", ...
    "X-averaged negative electrode SEI thickness [m]", ...
    "X-averaged negative electrode total SEI thickness [m]", ...
    "Negative SEI thickness [m]", ...
    "Negative total SEI thickness [m]"});

j_py_A = get_pybamm_var(sol, { ...
    "X-averaged negative SEI interfacial current density [A.m-2]", ...
    "X-averaged negative electrode SEI interfacial current density [A.m-2]", ...
    "X-averaged SEI interfacial current density [A.m-2]", ...
    "Negative SEI interfacial current density [A.m-2]", ...
    "Negative electrode SEI interfacial current density [A.m-2]"});

V_py = get_pybamm_var(sol, { ...
    "Terminal voltage [V]", ...
    "Voltage [V]"});

I_py = get_pybamm_var(sol, { ...
    "Current [A]", ...
    "Current variable [A]", ...
    "Cell current [A]"});

F = 96485.33212;
j_py = j_py_A ./ F;              % convert A/m2 to mol/m2/s
dLdt_py = gradient(L_py, t_py);

%% Numerical comparison
L_py_i = interp1(t_py, L_py, t_my, "linear", "extrap");
dLdt_py_i = interp1(t_py, dLdt_py, t_my, "linear", "extrap");
j_py_i = interp1(t_py, j_py, t_my, "linear", "extrap");

fprintf("\nInitial SEI thickness:\n");
fprintf("  My SPMe Ramadass       : %.6e m\n", L_my(1));
fprintf("  PyBaMM reaction limited: %.6e m\n", L_py(1));

fprintf("\nFinal values:\n");
fprintf("  My SPMe L_sei          = %.6e m\n", L_my(end));
fprintf("  PyBaMM L_sei           = %.6e m\n", L_py(end));
fprintf("  My SPMe dL/dt          = %.6e m/s\n", dLdt_my(end));
fprintf("  PyBaMM dL/dt           = %.6e m/s\n", dLdt_py(end));
fprintf("  My SPMe j_sei          = %.6e mol/m2/s\n", j_my(end));
fprintf("  PyBaMM j_sei/F         = %.6e mol/m2/s\n", j_py(end));

fprintf("\nRMSE:\n");
fprintf("  L_sei = %.6e m\n", sqrt(mean((L_my - L_py_i).^2)));
fprintf("  dLdt  = %.6e m/s\n", sqrt(mean((dLdt_my - dLdt_py_i).^2)));
fprintf("  j_sei = %.6e mol/m2/s\n", sqrt(mean((j_my - j_py_i).^2)));

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
title(sprintf("Terminal voltage comparison, %.2gC", C_rate));
legend("My SPMe Ramadass", "PyBaMM", "Location", "best");

figure;
plot(t_my/3600, I_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, I_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("Current [A]");
title(sprintf("Current comparison, %.2gC", C_rate));
legend("My SPMe Ramadass", "PyBaMM", "Location", "best");

figure;
plot(t_my/3600, L_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, L_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SEI thickness [m]");
title(sprintf("SEI thickness comparison, %.2gC", C_rate));
legend("My SPMe Ramadass", "PyBaMM reaction limited", "Location", "best");

figure;
plot(t_my/3600, dLdt_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, dLdt_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SEI growth rate dL/dt [m/s]");
title(sprintf("SEI growth rate comparison, %.2gC", C_rate));
legend("My SPMe Ramadass", "PyBaMM reaction limited", "Location", "best");

figure;
plot(t_my/3600, j_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, j_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("j_{SEI} [mol m^{-2} s^{-1}]");
title(sprintf("j_{SEI} comparison, %.2gC", C_rate));
legend("My SPMe Ramadass", "PyBaMM j_{SEI}/F", "Location", "best");

%% Local functions
function try_update_param(param, candidate_keys, val)
    for i = 1:numel(candidate_keys)
        key = candidate_keys{i};
        try
            param.update(py.dict(pyargs(key, val)), pyargs("check_already_exists", true));
            fprintf("Updated PyBaMM parameter: %s = %.6g\n", key, val);
            return;
        catch
        end
    end
    fprintf("Skipped PyBaMM parameter aliases: %s\n", strjoin(candidate_keys, " | "));
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
            fprintf("Failed variable: %s\n", name);
            fprintf("  %s\n", ME.message);
        end
    end

    error("Could not find requested PyBaMM variable. Tried aliases: %s", ...
        strjoin(cellstr(string(candidate_names)), " | "));
end

%% test_spme_safari_vs_pybamm_solvent_diffusion_limited.m
% Compare your SPMe.m Safari SEI model vs PyBaMM solvent-diffusion-limited SEI.
%
% Mapping:
%   Your model:       DegradationModel = "safari"
%   PyBaMM model:     options("SEI") = "solvent-diffusion limited"
%
% Note:
%   Your Safari branch solves dynamic EC diffusion through the SEI with a
%   reaction boundary condition. PyBaMM's "solvent-diffusion limited" option
%   is the reduced diffusion-limited form. Therefore, this test compares the
%   Safari diffusion-limited PyBaMM option against your more detailed Safari
%   implementation, not an equation-identical implementation.
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

initial_sei_m = 1.0e-6;     % force both models to same initial SEI (warm start at diffusion-limited regime for safari)
Nsei = 10;                  % must match SPMe.m default unless you modified it

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

%% Build initial SPMe state with matched initial SEI and EC-in-SEI profile
disable_thermal = true;
degradation_model = "safari";
electrolyte_model = "fdm";

p = load_SPMe_params(init_soc, Tamb_C, Nr, dt, Nxn, Nxs, Nxp, ...
    disable_thermal, degradation_model, electrolyte_model);

% Ensure initial state length matches SPMe.m Safari state order:
% [c_s_n; c_s_p; c_e; delta_sei; c_ec_sei]
p.Nsei = Nsei;
if ~isfield(p, "epsilon_sei") || isempty(p.epsilon_sei)
    p.epsilon_sei = 1.0;
end

calc = calculate_equilibrium_voltage_from_soc(init_soc);
[csn0, csp0] = init_cs(p, calc.Veq);

c_n0 = csn0 * ones(p.Nr - 1, 1);
c_p0 = csp0 * ones(p.Nr - 1, 1);
c_e0 = p.c_e * ones(p.Nx - 3, 1);

% Start EC concentration inside SEI at the outer-boundary value.
% c_ec_sei0 = p.epsilon_sei * p.c0_ec * ones(Nsei, 1);
c_outer = p.epsilon_sei * p.c0_ec;

% Inner surface ~ 0 for KL/D >> 1
c_s = 0;

z = (1:Nsei)'/(Nsei+1);

c_ec_sei0 = c_outer * z;


init_state = struct();
init_state.x = [c_n0; c_p0; c_e0; initial_sei_m; c_ec_sei0];


%% Run your SPMe.m Safari model
sim_my = sim_profile_LGM50("SPMe", init_soc, {"CC", C_rate}, ...
    "Nxn", Nxn, ...
    "Nxs", Nxs, ...
    "Nxp", Nxp, ...
    "Nr", Nr, ...
    "DeltaT", dt, ...
    "Tamb", Tamb_C, ...
    "DegradationModel", "safari", ...
    "ElectrolyteModel", "fdm", ...
    "UseElectrolyteActivity", false, ...
    "InitialState", init_state, ...
    "SPMeOutputs", {"delta_sei","delta_sei_dot","j_sei","j_side", ...
                    "eta_n_side","phi_se","eta_s","c_s_ec", ...
                    "safari_pe","safari_k","V","SOC"}, ...
    "SuppressSummary", true);

keep = sim_my.time <= T_end_s;

t_my = sim_my.time(keep);
L_my = sim_my.spme_outputs.delta_sei(keep);
dLdt_my = sim_my.spme_outputs.delta_sei_dot(keep);
j_my = sim_my.spme_outputs.j_sei(keep);   % mol/m2/s in your code
V_my = sim_my.voltage(keep);
I_my = sim_my.current(keep);

K_my = sim_my.spme_outputs.safari_k(keep);   % kinetic prefactor
KL_over_D = K_my .* L_my ./ p.D_ec;

%% Run PyBaMM from MATLAB Python interface
pybamm = py.importlib.import_module("pybamm");
np = py.importlib.import_module("numpy");

options = py.dict(pyargs( ...
    "SEI", "solvent-diffusion limited", ...
    "SEI film resistance", "distributed", ...
    "thermal", "isothermal"));

model = pybamm.lithium_ion.SPMe(options);
param = pybamm.ParameterValues(pybamm_param_set);

% Match key SEI parameters to your params_LGM50.m where possible.
% Your Safari branch uses:
%   p.D_ec, p.c0_ec, p.k0_sei, p.kappa_sei, p.M_sei, p.rho_sei
try_update_param(param, {"Initial SEI thickness [m]"}, initial_sei_m);
try_update_param(param, {"SEI solvent diffusivity [m2.s-1]", ...
                         "SEI solvent diffusivity [m^2.s-1]"}, p.D_ec);
try_update_param(param, {"Bulk solvent concentration [mol.m-3]", ...
                         "EC initial concentration in electrolyte [mol.m-3]", ...
                         "Solvent concentration [mol.m-3]"}, p.c0_ec);
try_update_param(param, {"SEI resistivity [Ohm.m]"}, 1 / p.kappa_sei);
try_update_param(param, {"SEI partial molar volume [m3.mol-1]"}, p.M_sei / (2 * p.rho_sei));

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
fprintf("  My SPMe Safari                  : %.6e m\n", L_my(1));
fprintf("  PyBaMM solvent-diffusion limited: %.6e m\n", L_py(1));

fprintf("\nFinal values:\n");
fprintf("  My SPMe L_sei                  = %.6e m\n", L_my(end));
fprintf("  PyBaMM L_sei                   = %.6e m\n", L_py(end));
fprintf("  My SPMe dL/dt                  = %.6e m/s\n", dLdt_my(end));
fprintf("  PyBaMM dL/dt                   = %.6e m/s\n", dLdt_py(end));
fprintf("  My SPMe j_sei                  = %.6e mol/m2/s\n", j_my(end));
fprintf("  PyBaMM j_sei/F                 = %.6e mol/m2/s\n", j_py(end));

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
legend("My SPMe Safari", "PyBaMM", "Location", "best");

figure;
plot(t_my/3600, I_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, I_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("Current [A]");
title(sprintf("Current comparison, %.2gC", C_rate));
legend("My SPMe Safari", "PyBaMM", "Location", "best");

figure;
plot(t_my/3600, L_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, L_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SEI thickness [m]");
title(sprintf("SEI thickness comparison, %.2gC", C_rate));
legend("My SPMe Safari", "PyBaMM solvent-diffusion limited", "Location", "best");

figure;
plot(t_my/3600, dLdt_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, dLdt_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SEI growth rate dL/dt [m/s]");
title(sprintf("SEI growth rate comparison, %.2gC", C_rate));
legend("My SPMe Safari", "PyBaMM solvent-diffusion limited", "Location", "best");

figure;
plot(t_my/3600, j_my, "LineWidth", 1.5); hold on;
plot(t_py/3600, j_py, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("j_{SEI} [mol m^{-2} s^{-1}]");
title(sprintf("j_{SEI} comparison, %.2gC", C_rate));
legend("My SPMe Safari", "PyBaMM j_{SEI}/F", "Location", "best");

figure;
plot(t_my/3600, KL_over_D, "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("KL/D [-]");
title("Diffusion-limited indicator (KL/D)");

yline(1, '--k', 'KL/D = 1');
yline(10, '--r', 'KL/D = 10');
yline(100, '--g', 'KL/D = 100');


dL_my_inc = L_my - L_my(1);
dL_py_inc = L_py - L_py(1);

figure;
plot(t_my/3600, dL_my_inc, "LineWidth", 1.5); hold on;
plot(t_py/3600, dL_py_inc, "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("\Delta SEI thickness [m]");
title("SEI growth increment comparison");
legend("My Safari", "PyBaMM diffusion-limited");

j_ratio = j_my ./ j_py_i;

figure;
plot(t_my/3600, j_ratio, "LineWidth", 1.5); hold on;
plot(t_my/3600, (KL_over_D ./ (1 + KL_over_D)), "--", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("Ratio");
title("j_{Safari} / j_{PyBaMM} vs theoretical limit");
legend("Observed", "KL/D / (1+KL/D)");

Rsei_my_voltage = L_my ./ p.kappa_P;      % what your current voltage model uses
Rsei_my_safari  = L_my ./ p.kappa_sei;    % what Safari physics would imply
Rsei_py         = L_py .* (1 / p.kappa_sei);

figure;
plot(t_my/3600, Rsei_my_voltage, "LineWidth", 1.5); hold on;
plot(t_my/3600, Rsei_my_safari, "--", "LineWidth", 1.5);
plot(t_py/3600, Rsei_py, ":", "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("SEI resistance [Ohm m^2]");
title("SEI film resistance comparison");
legend("My SPMe voltage uses L/kappa_P", ...
       "My Safari expected L/kappa_sei", ...
       "PyBaMM L * resistivity", ...
       "Location", "best");

j_diff_limit_my = -p.D_ec .* (p.epsilon_sei * p.c0_ec) ./ L_my;  % mol/m2/s
ratio_to_diff_limit = j_my ./ j_diff_limit_my;

figure;
plot(t_my/3600, ratio_to_diff_limit, "LineWidth", 1.5);
grid on;
xlabel("Time [h]");
ylabel("j_{Safari} / (-D c_0 / L)");
title("Safari flux relative to diffusion-limited maximum");
yline(1, "--");

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

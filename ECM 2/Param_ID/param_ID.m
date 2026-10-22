clc; clear; close all
addpath("..")
%% Set Model Order 
n_order = 2;

%% Load SPMe pulse test data
% load("spme_pulse.mat")
% train_duration = 1:100;
% if numel(train_duration) > 1
%     V_meas = spme.volt(train_duration);  % Terminal voltage (V)
%     I = spme.cur(train_duration);        % Current (A), positive = discharge
%     tdata = spme.time(train_duration);   % Time vector (s)
% else
%     V_meas = spme.volt;
%     I = spme.cur;
%     tdata = spme.time;
% end
load("dfn_-1C_25deg_0.1SOC.mat")
V_meas = dfn.volt;
I = dfn.cur;
tdata = dfn.time;

% Load OCV-SOC lookup data
load("..\..\param\OCVdata.mat", "ocv_p")
% Load Battery Params
p = {};
p.OneC = 5;


%% Initial States
% Estimate initial SOC by interpolating first voltage
% z0_calc = interp1(OCV_data, z_ocv, V_meas(1), 'linear', 'extrap');
% z0_calc = max(min(z0_calc, 1), 0);  % Clamp to [0,1]
ocv_p_inv = ocv_p;
ocv_p_inv(end) = ocv_p_inv(end) - V_meas(1);
z0_calc = real(roots(ocv_p_inv.'));
z0_calc = z0_calc(1);

% Initial state: [SOC; Vc]
x0 = zeros(n_order+1, 1);
% x0(1) = z0_calc;
x0(1) = dfn.soc(1);

% Define PSO bounds for [R0, R1, C1]
% lb = [0, 1e-3, 100,, ...]
lb = zeros(n_order*2+1, 1);
lb(2:2:end) = 1e-3;
lb(3:2:end) = 100;
% ub = [0.1, 0.07, 7000, ...]
ub = zeros(n_order*2+1, 1);
ub(1) = 0.1;
ub(2:2:end) = 0.07;
ub(3:2:end) = 7000;


% fprintf('Running Particle Swarm Optimization (PSO)...\n');
% options = optimoptions("particleswarm", ...
%     "Display", "iter", ...
%     "UseParallel", true, ...
%     'UseVectorized', false, ...
%     "SwarmSize", 200, ...
%     "FunctionTolerance", 1e-6);
% 
% % Run PSO
% params = particleswarm(@(p) errorFunction(p, tdata, I, V_meas, x0, OCV_data, z_ocv), ...
%     numel(lb), lb, ub, options);

fun = @(ecm_p) errorFunction(ecm_p, tdata, I, V_meas, x0, p, ocv_p, n_order);
nvars = numel(lb);
fprintf('Running Particle Swarm Optimization (PSO)...\n');
options = optimoptions("particleswarm","Display","iter", ...
            'UseParallel', true, ...
            'UseVectorized', false, ...
            "FunctionTolerance",1e-6,"SwarmSize",200);
ecm_p = particleswarm(fun,nvars,lb,ub,options);


% Display results
fprintf('\nIdentified Parameters:\n');
fprintf('R0 = %.6f Ohm\n', ecm_p(1));
for i = 1:n_order
    fprintf('R%d = %.6f Ohm\n', i, ecm_p(i*2)); 
    fprintf('C%d = %.2f F\n', i, ecm_p(i*2 + 1));
end
%% Save ECM Parameters
save("../../param/ECM_Params_order" + string(n_order) + ".mat", "ecm_p");

%% Solve ECM with optimized params
% u = @(t) interp1(spme.time, spme.cur, t, 'spline', 0);
% opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
% [tsol, xsol] = ode15s(@(t,x) ode_ecm(x, u(t), p, ecm_p, n_order), spme.time, x0);%, opts);
% V_sim = ecm_voltage(xsol, spme.cur, ecm_p, ocv_p, n_order);
% 
% % Compute and display RMSE
% rmse = sqrt(mean((spme.volt.' - V_sim).^2));
% fprintf('RMSE = %.6f V\n', rmse);
% 
% % Plot result
% figure;
% plot(spme.time, spme.volt, 'r-', 'LineWidth', 1.5); hold on;
% plot(tsol, V_sim, 'b--', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Voltage (V)');
% legend('Measured', 'Simulated', 'Location', 'southeast');
% title('RMSE = ' + string(rmse) + ' V');
% grid on;
% 
% savefig("ECM_fit_order" + string(n_order) + ".fig")

u = @(t) interp1(dfn.time, dfn.cur, t, 'spline', 0);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[tsol, xsol] = ode15s(@(t,x) ode_ecm(x, u(t), p, ecm_p, n_order), dfn.time, x0);%, opts);
V_sim = ecm_voltage(xsol.', dfn.cur, ecm_p, ocv_p, n_order);

% Compute and display RMSE
rmse = sqrt(mean((dfn.volt - V_sim).^2));
fprintf('RMSE = %.6f V\n', rmse);

% Plot result
figure;
plot(dfn.time, dfn.volt, 'r-', 'LineWidth', 1.5); hold on;
plot(tsol, V_sim, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Measured', 'Simulated', 'Location', 'southeast');
title('RMSE = ' + string(rmse) + ' V');
grid on;

savefig("ECM_fit_order" + string(n_order) + ".fig")


%% 
rmpath("..")
% === SUBFUNCTIONS ===
function err = errorFunction(ecm_p, tdata, I, V_meas, x0, p, ocv_p, n_order)
    u = @(t) interp1(tdata, I, t, 'spline', 0);
%     opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [~, xsol] = ode15s(@(t,x) ode_ecm(x, u(t), p, ecm_p, n_order), tdata, x0); %, opts);
    V_sim = ecm_voltage(xsol.', I, ecm_p, ocv_p, n_order);

    err = sqrt(mean((V_meas - V_sim).^2));
end


clear;
clc;
% close all;
addpath("DFN")
addpath("param")

disp('Doyle-Fuller-Newman Model (DFN)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation HyperParameters
verbosity = 0;
visualize_ocv_curve = true;
save_data = true;
disable_thermal = true; disable_sei = true;

%% Loop to collect data
CRate = 1/100;
t_amb = 25;

%% DFN Params
delta_t = 1;
load_init = false;

p = load_DFN_params(0, t_amb, load_init, delta_t, disable_thermal, disable_sei);

%% Input Discharge Current
I = p.OneC * CRate; 

% Estimate Simulation length to initialize Vectors
NT = floor(3800/abs(CRate))*5;

%% Initialize States & Preallocation
x0 = p.x0;
z0 = p.z0;
y0 = p.y0;

x = zeros(length(x0), NT);
x(:, 1) = x0;
z = zeros(length(z0), NT);
z(:, 1) = z0;
y = zeros(length(y0), NT);
y(:, 1) = y0;
Nn = p.Nxn - 1; Np = p.Nxp - 1; Nx = p.Nx - 3;
%% Integrate
if verbosity > 0 
    disp('Simulating DFN Model...')
end

%% Charge Battery to Max Voltage
Volt = 4.183;
k = 0;

while Volt < p.volt_max
    T = x0(end);
    % Current
    Cur_vec = [-I -I -I];

    % Step
    [x(:, 1), z(:, 1), y(:, 1), stats] = step_dfn(x0, z0, Cur_vec, p, verbosity); 
    Volt = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    z0 = z(:, 1);
    if verbosity > 1
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
        k*p.delta_t,-I/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
end
%% Hold for 1.5hr
V_last = Volt;
for k = 1:5400
    T = x0(end);
    I_charge = max([((Volt - p.volt_max) *  10 + ...
                    (Volt - V_last) *  0.01);
                   -I]);
    V_last = Volt;
    % Current
    Cur_vec = [Cur_vec(2:end), I_charge];

    % Step
    [x(:, 1), z(:, 1), y(:, 1), stats] = step_dfn(x0, z0, Cur_vec, p, verbosity); 
    Volt = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    z0 = z(:, 1);

    if verbosity > 1
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
        k*p.delta_t,I_charge/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
end
%% Rest for 30min
for k = 1:1800
    T = x0(end);
    % Current
    Cur_vec = [0 0 0];

    % Step
    [x(:, 1), z(:, 1), y(:, 1), stats] = step_dfn(x0, z0, Cur_vec, p, verbosity); 
    Volt = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    z0 = z(:, 1);

    if verbosity > 1
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
        k*p.delta_t,0,T-273.15,SOC,Volt,stats.iters);
    end
end

%% Discharge Battery
k = 0;
while SOC > 1e-3
Cur_vec = [0, 0, 0];
% for k = 1:NT
    k = k+1;
    T = x(end, k);
    % Current
    I_discharge = min([((Volt - p.volt_min) *  10 + ...
                    (Volt - V_last) *  0.01);
                   I]);
    V_last = Volt;
    % Current
    Cur_vec = [Cur_vec(2:end), I_discharge];

    % Step
    [x(:, k+1), z(:, k+1), y(:, k+1), stats] = step_dfn(x(:,k),z(:,k), Cur_vec, p, verbosity);
    Volt = y(1, k+1);
    SOC = y(2, k+1);
    c_ex = y(3*Nn + 2*Np + 3 : 3*Nn + 2*Np + Nx + 6, k+1);

    if verbosity > 1
        fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
            k*p.delta_t,I_discharge/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
    t_last = k;
end

%% Parse out states
[Volt, SOC, ~, ~, ~, ~, ~, ~, ~, ~] = parse_y(y(:, 1:k), t_last, p);
ocv_discharge = Volt.';
z_discharge = SOC.';

%% Reinitialize States for charging
x0 = x(:, k);
y0 = y(:, k);
z0 = z(:, k);
x = zeros(length(x0), NT);
x(:, 1) = x0;
z = zeros(length(z0), NT);
z(:, 1) = z0;
y = zeros(length(y0), NT);
y(:, 1) = y0;

%% Rest for 30min

for k = 1:1800
    T = x0(end);
    % Current
    Cur_vec = [0 0 0];

    % Step
    [x(:, 1), z(:, 1), y(:, 1), stats] = step_dfn(x0, z0, Cur_vec, p, verbosity); 
    Volt = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    z0 = z(:, 1);
    if verbosity > 1
    fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
        k*p.delta_t,0,T-273.15,SOC,Volt,stats.iters);
    end
end

%% Charge Battery
k = 0;
while SOC < 1 - 1e-3
    k = k+1;
    T = x(end, k);

    % Current
    I_charge = max([((Volt - p.volt_max) *  10 + ...
                    (Volt - V_last) *  0.01);
                   -I]);
    V_last = Volt;
    % Current
    Cur_vec = [Cur_vec(2:end), I_charge];

    % Step
    [x(:, k+1), z(:, k+1), y(:, k+1), stats] = step_dfn(x(:,k),z(:,k), Cur_vec, p, verbosity); 
    Volt = y(1, k+1);
    SOC = y(2, k+1);
    c_ex = y(3*Nn + 2*Np + 3 : 3*Nn + 2*Np + Nx + 6, k+1);
    
    if verbosity > 1
        fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
            k*p.delta_t,I_charge/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
    t_last = k;
end

%% Parse out states
[Volt, SOC, ~, ~, ~, ~, ~, ~, ~, ~] = parse_y(y(:, 1:k), t_last, p);
ocv_charge = Volt.';
z_charge = SOC.';

%% Fit OCV curve
z = linspace(0, 1, 1000);

ocv_discharge_interp = interp1(z_discharge, ocv_discharge, z, 'nearest', 'extrap');
ocv_charge_interp = interp1(z_charge, ocv_charge, z, 'nearest', 'extrap');
ocv_avg = 0.5 * (ocv_discharge_interp + ocv_charge_interp);

n = 9;
ocv_p = polyfit(z, ocv_avg, n);
ocv_fit = polyval(ocv_p, z);

if visualize_ocv_curve
    figure
    clf
    hold on
    plot(z_charge, ocv_charge, z_discharge, ocv_discharge, z, ocv_fit)
    xlabel('SOC')
    ylabel('OCV')
    legend("charge", "discharge", "Poly fit")
end
%% Save Data
beep

save("ECM\OCVdata.mat", "ocv_discharge", "z_discharge", "ocv_charge", "z_charge", "ocv_p", "n");
%% Restore Path
rmpath("DFN")
rmpath("param")
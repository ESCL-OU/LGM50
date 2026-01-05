clear;
clc;
close all;
import casadi.*
addpath("ECM")
addpath("param")

disp('Equivalent Circuit Model (ECM)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation HyperParameters
verbosity = 3; delta_t = 1;
save_data = true;
n_order = 2;

%% Loop to collect data
CRates = -1.5;
t_ambs = 25;
init_SOCs = 0.1;

for CRate = CRates
    for t_amb = t_ambs
        for init_SOC = init_SOCs
clc;

%% Input (dis)charge Current && Load ECM
dfn_dir = 'results\DFN\';
if CRate < 0
    dfn_file = "dfn_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOCs) + "SOC.mat";
else
    dfn_file = "dfn_" + string(CRate) + "C_" + string(t_amb) + "deg.mat";
end
if exist(dfn_dir + dfn_file, 'file')
    load(dfn_dir + dfn_file)
    t = dfn.time;
    I = dfn.cur;
    z0 = dfn.soc(1);
    [p, ecm_p, ocv_p] = load_ECM_params(z0, n_order);
    ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');
else
    if CRate == 0
        t = 0:delta_t:(100);
        I = zeros(size(t)); 
        z0 = init_SOC;
        [p, ecm_p, ocv_p] = load_ECM_params(z0, n_order);
        ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');
    else
        t = 0:delta_t:(floor(3600/abs(CRate))); 
        if CRate < 0 % Charging
            z0 = init_SOC;
        else % Discharging
            z0 = 1;
        end
        [p, ecm_p, ocv_p] = load_ECM_params(z0, n_order);
        ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');
        I = p.OneC*ones(size(t)) * CRate; 
    end
end
NT = length(t);

%% Initialize States & Preallocation
x0 = p.x0;
y0 = p.y0;

x = zeros(length(x0), NT);
x(:, 1) = x0;
V = zeros(length(y0), NT);
V(:, 1) = y0;

%% Integrate
if verbosity > 0
    disp('Simulating ECM ' + string(n_order) + "th order Model...")
end

for k = 1:(NT-1)
    [x(:, k+1), V(:, k+1)] = ecm_ss.step(x(:, k), I(k));
    SOC = x(1, k+1);
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,I(k)/p.OneC,p.T_amb-273.15, x(1, k+1),V(1, k+1));
    end
end

%% Save Results
ecm.time = t;
ecm.cur = I;
ecm.x = x;
ecm.volt = V;
if CRate < 0
    save("results\ECM\ecm^"+ string(n_order) + "_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOC) + "SOC.mat", ...
        'ecm');
else
    save("results\ECM\ecm^"+ string(n_order) + "_" + string(CRate) + "C_" + string(t_amb) + "deg.mat", ...
        'ecm');
end
        end
    end
end
beep

%% Restore Path
rmpath("ECM")
rmpath("param")
clear;
clc;
close all;
addpath("DFN")
addpath("param")

disp('Doyle-Fuller-Newman Model (DFN)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation HyperParameters
verbosity = 1;
save_data = true;
disable_thermal = true; disable_sei = true;

%% Loop to collect data
CRates = -0.1:-0.1:-1.5;
t_ambs = [25];
init_SOCs = [0.3];
target_SOC = 0.9;

for CRate = CRates
    for t_amb = t_ambs
        for init_SOC = init_SOCs
clc;

%% DFN Params
delta_t = 1;
if CRate < 0
    load_init = true;
else
    load_init = false;
end
p = load_DFN_params(init_SOC, t_amb, load_init, delta_t, disable_thermal, disable_sei);
if CRate < 0
    target_V = p.volt_max;
else
    target_V = p.volt_min;
end
%% Constant Current
NT = floor(3600/abs(CRate)) * 2;
t = 0:p.delta_t:NT; 
I = p.OneC * ones(size(t)) * CRate;

%% Constant Voltage Control
kp = 27;
kd = 2;
ki = 1.2;

epsilon = 0.3 * sign(CRate);
eps_decay = 0.94;

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
    disp('init. SOC: ' + string(init_SOC) + " | target SOC: " + string(target_SOC))
    disp('C-Rate: ' + string(CRate))
end

Volt = y(1, 1);
V_last = Volt;

%% Constant Current
SOC = init_SOC;
k = 1;
while SOC <= target_SOC
    T = x(end, k);
    % Current
    I(k:k+1) = CRate * p.OneC;
    V_last = Volt;
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
    % Step
    [x(:, k+1), z(:, k+1), y(:, k+1), stats] = step_dfn(x(:,k),z(:,k), Cur_vec, p, verbosity); 
    Volt = y(1, k+1);
    SOC = y(2, k+1);

    if verbosity > 1
        fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
            t(k),I(k+1)/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
    t_last = k;
    if(Volt < p.volt_min)
        break;
    elseif(Volt > p.volt_max)
        break;
    end
    k = k + 1;
end
time_cc_end = t_last;

volt_last =  y(1, k);
int_error = epsilon;
%% Constant Voltage
while SOC <= target_SOC
    T = x(end, k);
    % Current
    error = Volt - target_V + epsilon;
    d_error = (Volt - volt_last)  * p.delta_t;
    int_error = int_error + error * p.delta_t;
    volt_last = Volt;
    if CRate < 0
        if SOC >= target_SOC
            I(k+1) = 0;
        else
            I(k+1) = min([max([error * kp + ...
                        d_error * kd + ...
                        int_error * ki;
                        CRate*p.OneC]), 0]);   
        end
    else
        if SOC <= target_SOC
            I(k+1) = 0;
        else
            I(k+1) = max([min([error * kp + ...
                        d_error * kd + ...
                        int_error * ki;
                        CRate*p.OneC]), 0]);
        end
    end
    epsilon = epsilon * eps_decay;
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
    % Step
    [x(:, k+1), z(:, k+1), y(:, k+1), stats] = step_dfn(x(:,k),z(:,k), Cur_vec, p, verbosity); 
    Volt = y(1, k+1);
    SOC = y(2, k+1);
    
    if verbosity > 1
        fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
            t(k),I(k+1)/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
    t_last = k;
    k = k + 1;
end
%% Parse out states
[c_s_n, c_s_p, c_e, T] = parse_x(x, t_last, p);
[phi_s_n, phi_s_p, i_en, i_ep, phi_e, jn, jp] = parse_z(z, t_last, p);
[Volt, SOC, eta_pl, c_ss_n, c_ss_p, c_avg_n, c_avg_p, c_ex, eta_n, eta_p, Unref] = parse_y(y, t_last, p);

%% Save results
dfn.time = t(1:t_last);
dfn.cur = I(1:t_last);
dfn.c_s_n = c_s_n;
dfn.c_s_p = c_s_p;
dfn.c_e = c_e;
dfn.T = T;

dfn.phi_s_n = phi_s_n;
dfn.phi_s_p = phi_s_p;
dfn.i_en = i_en;
dfn.i_ep = i_ep;
dfn.phi_e = phi_e;
dfn.jn = jn;
dfn.jp = jp;

dfn.volt = Volt;
dfn.soc = SOC;
dfn.eta_pl = eta_pl;
dfn.c_ss_n = c_ss_n;
dfn.c_ss_p = c_ss_p;
dfn.c_avg_n = c_avg_n;
dfn.c_avg_p = c_avg_p;
dfn.c_ex = c_ex;
dfn.eta_n = eta_n;
dfn.eta_p = eta_p;
dfn.Unref = Unref;

%%
if  save_data
    if CRate < 0
    save("results\DFN_CCCV\dfn_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOC) + "SOC0_" + string(target_SOC) + "SOC.mat", ...
        'dfn');
    else
        save("results\DFN_CCCV\dfn_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(target_SOC) + "SOC.mat", ...
        'dfn');
    end

end
        end
    end
end
beep

%% Restore Path
rmpath("DFN")
rmpath("param")
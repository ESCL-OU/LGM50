clear;
clc;
% close all;
addpath("DFN")
addpath("param")

disp('Doyle-Fuller-Newman Model (DFN)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation HyperParameters
verbosity = 3;
save_data = true;
disable_thermal = true; disable_sei = true;
% For debugging purposes
collect_potentials = true;

%% Loop to collect data
CRates = -0.1:-0.1:-1.5;
t_ambs = 0:40;
init_SOCs = [0.1, 0.15, 0.2, 0.3, 0.4];

CRates = -1;
t_ambs = [3, 25];
init_SOCs = 0.1;

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

%% Input (dis)charge Current
V0 = p.y0(1);
if CRate == 0
    t = 0:p.delta_t:(100);
    I = zeros(size(t)); 
else
    t = 0:p.delta_t:(floor(3600/abs(CRate))); 
    I = p.OneC*ones(size(t)) * CRate; 
end
NT = length(t);

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

% Collect Potentials
if collect_potentials
    eta_n = zeros(Nn, NT);
    eta_p = zeros(Np, NT);
    Unref = zeros(Nn, NT);
    Upref = zeros(Np, NT);
    IR_n = zeros(Nn, NT);
    IR_p = zeros(Np, NT);
    phi_e_n = zeros(Nn, NT);
    phi_e_p = zeros(Np, NT);
    V_noVCE = zeros(1, NT);
    V_electrolyteCond = zeros(1, NT);
    V_electrolytePolar = zeros(1, NT);
    [eta_n(:, 1), eta_p(:, 1), Unref(:, 1), Upref(:, 1), ...
     IR_n(:, 1), IR_p(:, 1), phi_e_n(:, 1), phi_e_p(:, 1), ...
     V_noVCE(1), V_electrolyteCond(1), V_electrolytePolar(1)] ...
                = potentials(x(:, 1), z(:, 1), 0, p);
end


%% Integrate
if verbosity > 0 
    disp('Simulating DFN Model...')
end


for k = 1:(NT - 1)
    T = x(end, k);
    % Current
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end

    % Step
    [x(:, k+1), z(:, k+1), y(:, k+1), stats] = step_dfn(x(:,k),z(:,k), Cur_vec, p, verbosity); 
    Volt = y(1, k+1);
    SOC = y(2, k+1);
    c_ex = y(3*Nn + 2*Np + 3 : 3*Nn + 2*Np + Nx + 6, k+1);
    
    if collect_potentials
        [eta_n(:, k+1), eta_p(:, k+1), Unref(:, k+1), Upref(:, k+1), ...
         IR_n(:, k+1), IR_p(:, k+1), phi_e_n(:, k+1), phi_e_p(:, k+1), ...
         V_noVCE(k+1), V_electrolyteCond(k+1), V_electrolytePolar(k+1)] ...
                    = potentials(x(:, k+1), z(:, k+1), I(k), p);
    end

    if verbosity > 1
        fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
            t(k),I(k+1)/p.OneC,T-273.15,SOC,Volt,stats.iters);
    end
    if ~isa(CRate,'string')
        if(Volt < p.volt_min)
            fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);
    %         beep;
            t_last = k; %YK modified(added line)
            break; %YK modified(activated line)
        elseif(Volt > p.volt_max)
            fprintf(1,'Max Voltage of %1.1fV exceeded\n',p.volt_max);
    %         beep;
            t_last = k; %YK modified(added line)
            break;
        elseif(any(c_ex < 1))
            fprintf(1,'c_e depleted below 1 mol/m^3\n');
    %         beep;
            t_last = k; %YK modified(added line)
            break; %YK modified(activated line)
        end
    end
    t_last = k;
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
dfn.c_en = c_e(1:p.Nxn, :);
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

if collect_potentials
    dfn.eta_n = eta_n(:, 1:t_last);
    dfn.eta_p = eta_p(:, 1:t_last);
    dfn.Unref = Unref(:, 1:t_last);
    dfn.Upref = Upref(:, 1:t_last);
    dfn.IR_n = IR_n(:, 1:t_last);
    dfn.IR_p = IR_p(:, 1:t_last);
    dfn.phi_e_n = phi_e_n(:, 1:t_last);
    dfn.phi_e_p = phi_e_p(:, 1:t_last);
    dfn.V_noVCE = V_noVCE(1:t_last);
    dfn.V_electrolyteCond = V_electrolyteCond(1:t_last);
    dfn.V_electrolytePolar = V_electrolytePolar(1:t_last);
end


if CRate < 0
    save("results\DFN\dfn_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOC) + "SOC.mat", ...
        'dfn');
else
    save("results\DFN\dfn_" + string(CRate) + "C_" + string(t_amb) + "deg.mat", ...
    'dfn');
end
        end
    end
end
beep

%% Restore Path
rmpath("DFN")
rmpath("param")
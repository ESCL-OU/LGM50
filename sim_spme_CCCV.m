%% Single Particle Model w/ Electrolyte & Temperature simulated with CasADi
%   Published December 18, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/
%
%   Code based on publications
%   Battery State Estimation for a Single Particle Model with Electrolyte Dynamics 
%   S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic 
%   IEEE Transactions on Control System Technology, to appear 
%   DOI: 10.1109/TCST.2016.2571663
%
%   Optimal Charging of Batteries via a Single Particle Model with Electrolyte and Thermal Dynamics 
%   H. Perez, X. Hu, S. J. Moura 
%   2016 American Control Conference
%   DOI: 10.1109/ACC.2016.7525538
%
%   Modified Feb 7, 2025 by Yuichi Kajiura (line specified as "YK modified")
%   Modified Feb 24, 2025 by Wonoo Choo to be modelled using CasADi
% clear;
clc;
% close all;
import casadi.*
addpath("SPMe")
addpath("param")

disp('Single Particle Model w/ Electrolyte (SPMe)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation HyperParameters
verbosity = 0;
save_data = true;
disable_thermal = true; disable_sei = true;
% For debugging purposes
collect_potentials = true;

%% Loop to collect data
CRates = -0.1:-0.1:-1.5;
t_ambs = [25];
init_SOCs = [0.3];
target_SOC = 0.9;

for CRate = CRates
    for t_amb = t_ambs
        for init_SOC = init_SOCs
clc;
% clearvars -except CRates t_ambs CRate t_amb save_data verbosity;

if CRate < 0
    load_init = true;
else
    load_init = false;
end

%% Load SPMe Params
Nr = 30;
Nxn = 70; Nxs = 35; Nxp = 70;
delta_t = 1;
p = load_SPMe_params(init_SOC, t_amb, load_init, Nr, delta_t, Nxn, Nxs, Nxp, disable_thermal, disable_sei);
x0 = p.x0;
spme_ss = SPMe(p, 'idas');

%%
%%%%%%%%%%%%%%% MANUAL INPUT WITH C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
dfn_dir = 'results\DFN_CCCV\';
if CRate < 0
    dfn_file = "dfn_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOC) + "SOC0_" + string(target_SOC) + "SOC.mat";
else
    dfn_file = "dfn_" + string(CRate) + "C_" + string(t_amb) + "deg.mat";
end
if exist(dfn_dir + dfn_file, 'file')
    load(dfn_dir + dfn_file)
    t = dfn.time;
    I = dfn.cur;
else
    if CRate == 0
        t = 0:p.delta_t:(100);
        I = zeros(size(t)); 
    else
        t = 0:p.delta_t:(floor(3600/abs(CRate))); 
        I = p.OneC*ones(size(t)) * CRate; 
    end
end
NT = length(t);

%% Initialize States & Preallocation
x = zeros(length(x0), NT);
y = zeros(length(spme_ss.y), NT);
c_ex = zeros(p.Nx+1, NT);
c_ss_n = zeros(1, NT);
c_ss_p = zeros(1, NT);
x(:, 1) = x0;
y0 = spme_ss.h('x', x(:, 1), 'p', 0);
y(:, 1) = full(y0.y);
Fss = spme_ss.f_ss('x', x(:, 1), 'u', 0);
c_ss_n(1) = full(evalf(Fss.c_ss_n));
c_ss_p(1) = full(evalf(Fss.c_ss_p));
Fcex = spme_ss.f_cex('x', x(:, 1));
c_ex(:, 1) = full(evalf(Fcex.c_ex));

% Collect Potentials
if collect_potentials
    eta_n = zeros(1, NT);
    eta_p = zeros(1, NT);
    Unref = zeros(1, NT);
    Upref = zeros(1, NT);
    V_noVCE = zeros(1, NT);
    V_electrolyteCond = zeros(1, NT);
    V_electrolytePolar = zeros(1, NT);
    [eta_n(1), eta_p(1), Unref(1), Upref(1), V_noVCE(1), V_electrolyteCond(1), V_electrolytePolar(1)] = ...
                    spme_ss.potentials(x(:, 1), 0);
end

%% Simulate SPMe Plant
tic;
disp('Simulating SPMe Plant...');
for k = 1:NT-1
    [x(:, k+1), y(:, k+1), c_ss_n(k+1), c_ss_p(k+1), c_ex(:, k+1)] = spme_ss.step(x(:, k), I(k));
    V = y(1, k+1);
    SOC = y(2, k+1);

    if collect_potentials
        [eta_n(k+1), eta_p(k+1), Unref(k+1), Upref(k+1), V_noVCE(k+1), V_electrolyteCond(k+1), V_electrolytePolar(k+1)] = ...
                    spme_ss.potentials(x(:, k+1), I(k));
    end
    
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,I(k)/p.OneC,p.T_amb-273.15,SOC,V);
    end
end
t_last = NT;

% Parse States
c_s_n = x(1:(p.Nr-1), 1:t_last);
c_s_p = x(p.Nr : 2*(p.Nr-1), 1:t_last);
c_e = x(2*p.Nr-1 : 2*p.Nr-1+p.Nx-4, 1:t_last);

V = y(1, 1:t_last);
SOC = y(2, 1:t_last);
eta_pl = y(3, 1:t_last);

spme.time = t;
spme.cur = I;
spme.c_s_n = c_s_n;
spme.c_s_p = c_s_p;
spme.c_e = c_e;
spme.c_en = c_e(1:p.Nxn, :);
spme.c_ex = c_ex;
spme.c_ss_n = c_ss_n;
spme.c_ss_p = c_ss_p;
spme.volt = V;
spme.soc = SOC;
spme.eta_pl = eta_pl;

if collect_potentials
    spme.eta_n = eta_n(1:t_last);
    spme.eta_p = eta_p(1:t_last);
    spme.Unref = Unref(1:t_last);
    spme.Upref = Upref(1:t_last);
    spme.V_noVCE = V_noVCE;
    spme.V_electrolyteCond = V_electrolyteCond(1:t_last);
    spme.V_electrolytePolar = V_electrolytePolar(1:t_last);
end

%%
if CRate < 0
    save("results\SPMe_CCCV\spme_" + string(CRate) + "C_" + string(t_amb) + "deg_" + string(init_SOC) + "SOC0_" + string(target_SOC) + "SOC.mat", ...
        'spme');
else
    save("results\SPMe_CCCV\spme_" + string(CRate) + "C_" + string(t_amb) + "deg" + string(target_SOC) + "SOC.mat", ...
        'spme');
end
        end
    end
end
beep

%% Restore Path
rmpath("SPMe")
rmpath("param")
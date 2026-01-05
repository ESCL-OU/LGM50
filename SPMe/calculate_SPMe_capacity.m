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

clear;
clc;
% close all;
import casadi.*
verbosity = 0;
save_data = true;

disp('Single Particle Model w/ Electrolyte (SPMe)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Simulation Param
CRate = 0.01;
t_amb = 25;

%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LGM50

%% Set ambient temp and if use thermal/sei dynamics (YK modified/added)
% p.T_amb = 25 + 273.15; 
p.T_amb = t_amb + 273.15;
p.disable_thermal = true;
p.disable_sei = true;
%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;

p.c_n_soc0 = cn_low;
p.c_p_soc0 = cp_low;
p.c_n_soc100 = cn_high;
p.c_p_soc100 = cp_high;

%%
%%%%%%%%%%%%%%% MANUAL INPUT WITH C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
I = p.OneC * CRate; 

% Estimate Simulation length to initialize Vectors
NT = floor(3800/abs(CRate));

%% Preallocation & Initial Conditions

%%% Finite difference for spherical particle
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
r_vec = (0:p.delta_r_n:1)';
r_vecx = r_vec(2:end-1);

% Finite difference points along x-coordinate
p.Nx = p.Nxn+p.Nxs+p.Nxp;
Nx = p.Nx - 3;
x_vec_spme = linspace(0,1,Nx+4);


p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

% Output Discretization params
disp('Discretization Params:');
fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
fprintf(1,'No. of FDM nodes in Single Particles : %1.0f\n',p.Nr);
fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
disp(' ');


%%% INITIAL CONDITIONS
% Solid concentration
csn0 = p.csn0;
csp0 = p.csp0;
c_n0 = csn0 * ones(p.Nr-1,1);
c_p0 = csp0 * ones(p.Nr-1,1);

% Electrolyte concentration
ce0 = p.c_e;
c_e0 = ce0 * ones(Nx,1);

% Temperature
T0 = p.T_amb;

% SEI layer
disp('Initial Conditions:');
fprintf(1,'Voltage : %1.3f V\n',p.volt_max);
fprintf(1,'Normalized Solid Concentration in Anode | Cathode : %1.2f | %1.2f\n',csn0/p.c_s_n_max,csp0/p.c_s_p_max);
fprintf(1,'Electrolyte Concentration : %2.3f kmol/m^3\n',ce0/1e3);
fprintf(1,'Temperature : %3.2f Celsius\n',T0-273.15);
disp(' ');

%% Generate Constant System Matrices

% Electrolyte concentration matrices
[M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce;

%% Simulate SPMe Plant
tic;
% Initial Conditions
x0 = [c_n0; c_p0; c_e0];
p.x0 = x0;

disp('Simulating SPMe Plant...');
spme_ss = SPMe(p, 'idas');


x = zeros(length(x0), NT);
y = zeros(length(spme_ss.y), NT);
c_ex = zeros(p.Nx+1, NT);
c_ss_n = zeros(1, NT);
c_ss_p = zeros(1, NT);
eta_n = zeros(1, NT);
x(:, 1) = x0;

%% Charge Battery to Max Voltage
V = 4.183;
k = 0;
while V < p.volt_max
    Fx = spme_ss.f_sim_etan('x', x0, 'p', -I);
    x(:, 1) = full(evalf(Fx.xk_1));
    y(:, 1) = full(evalf(Fx.y));
    V = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,-I/p.OneC,p.T_amb-273.15,SOC,V);
    end
end
%% Hold for 1.5hr
V_last = V;
for k = 1:5400
    I_charge = max([(V - p.volt_max) *  10+ ...
             (V - V_last) *  0.01;
             -I]);
    V_last = V;
    Fx = spme_ss.f_sim_etan('x', x0, 'p', I_charge);
    x(:, 1) = full(evalf(Fx.xk_1));
    y(:, 1) = full(evalf(Fx.y));
    V = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,I_charge/p.OneC,p.T_amb-273.15,SOC,V);
    end
end

%% Rest for 30 min
for k = 1:1800
    Fx = spme_ss.f_sim_etan('x', x0, 'p', 0);
    x(:, 1) = full(evalf(Fx.xk_1));
    y(:, 1) = full(evalf(Fx.y));
    V = y(1, 1);
    SOC = y(2, 1);
    x0 = x(:, 1);
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,0/p.OneC,p.T_amb-273.15,SOC,V);
    end
end
%% Discharge Battery
k = 0;
while V > p.volt_min
    k = k+1;
    Fx = spme_ss.f_sim_etan('x', x(:,k), 'p', I);
    x(:, k+1) = full(evalf(Fx.xk_1));
    y(:, k) = full(evalf(Fx.y));
    eta_n(k) = full(Fx.eta_n);
    Fss = spme_ss.f_ss('x', x(:, k), 'u', I);
    c_ss_n(k) = full(evalf(Fss.c_ss_n));
    c_ss_p(k) = full(evalf(Fss.c_ss_p));
    Fcex = spme_ss.f_cex('x', x(:, k));
    c_ex(:, k) = full(evalf(Fcex.c_ex));
    V = y(1, k);
    SOC = y(2, k);
    if verbosity > 1
       fprintf(1,'Time : %3.2f sec | C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV \n',...
           k*p.delta_t,I/p.OneC,p.T_amb-273.15,SOC,V);
    end
end
fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);

% Parse States
c_s_n = x(1:(p.Nr-1), 1:k);
c_s_p = x(p.Nr : 2*(p.Nr-1), 1:k);
c_e = x(2*p.Nr-1 : 2*p.Nr-1+p.Nx-4, 1:k);

V = y(1, 1:k);
SOC = y(2, 1:k);
eta_pl = y(3, 1:k);

spme.time = 0:p.delta_t:(k-1)*p.delta_t;
spme.cur = I * ones(size(spme.time));
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
spme.eta_n = eta_n;
%%
save("results\SPMe\spme_" + string(CRate) + "C_" + string(t_amb) + "deg.mat", ...
    'spme');

%% Compute Capacity
capacity = I * k * p.delta_t / 3600
%%
beep
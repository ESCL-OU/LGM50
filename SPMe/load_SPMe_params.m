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

function p = load_SPMe_params(init_SOC, t_amb, Nr, delta_t, Nxn, Nxs, Nxp, disable_thermal, degradation_model, electrolyte_model)
%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LGM50

if ~isscalar(init_SOC) || init_SOC < 0 || init_SOC > 1
    error('init_SOC must be a scalar in [0, 1].');
end

%% Set ambient temp and thermal dynamics switch
p.T_amb = t_amb + 273.15;
p.disable_thermal = disable_thermal;
p.delta_t = delta_t;
if nargin < 9 || isempty(degradation_model)
    degradation_model = 'none';
end
p.degradation_model = lower(char(degradation_model));
if ~ismember(p.degradation_model, {'none', 'ramadass', 'yang'})
    error('degradation_model must be one of: none, ramadass, yang.');
end
if nargin < 10 || isempty(electrolyte_model)
    electrolyte_model = 'fdm';
end
p.electrolyte_model = lower(char(electrolyte_model));
if ~ismember(p.electrolyte_model, {'fdm', 'fvm'})
    error('electrolyte_model must be one of: fdm, fvm.');
end

%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);

p.c_n_soc0 = cn_low;
p.c_p_soc0 = cp_low;
p.c_n_soc100 = cn_high;
p.c_p_soc100 = cp_high;

%% Preallocation & Initial Conditions
%%% Finite difference for spherical particle
p.Nr = Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;

% Finite difference points along x-coordinate
p.Nxn = Nxn;
p.Nxs = Nxs;
p.Nxp = Nxp;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

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
repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_dir);
calc = calculate_equilibrium_voltage_from_soc(init_SOC);
[csn0, csp0] = init_cs(p, calc.Veq);
p.csn0 = csn0;
p.csp0 = csp0;
ce0 = p.c_e;

c_n0 = csn0 * ones(p.Nr-1,1);
c_p0 = csp0 * ones(p.Nr-1,1);
c_e0 = ce0 * ones(p.Nx-3, 1);
V0 = calc.Veq;

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;

% SEI layer
delta_sei0 = 0;

disp('Initial Conditions:');
fprintf(1,'Voltage : %1.3f V\n',V0);
fprintf(1,'Normalized Solid Concentration in Anode | Cathode : %1.2f | %1.2f\n',csn0/p.c_s_n_max,csp0/p.c_s_p_max);
fprintf(1,'Electrolyte Concentration : %2.3f kmol/m^3\n',ce0/1e3);
fprintf(1,'Temperature in Roll | Can : %3.2f K | %3.2f K \n',T10,T20);
fprintf(1,'SEI Layer in Anode : %2f um \n',delta_sei0*1e6);
disp(' ');

%% Initial Condition
% Initial Conditions
x0 = [c_n0; c_p0; c_e0];

if ~p.disable_thermal
    x0 = [x0; T10; T20];
end
if ~strcmpi(char(p.degradation_model), 'none')
    x0 = [x0; delta_sei0];
end
p.x0 = x0;

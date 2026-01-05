%% Params for Electrochemical Model
%   Created May 15, 2025 by Wonoo Choo
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell

p.OneC = 5;
%% Geometric Params
% Thickness of each layer
p.L_n = 8.52e-5;     % Thickness of negative electrode [m]
p.L_s = 1.2e-5;     % Thickness of separator [m]
p.L_p = 7.56e-5;     % Thickness of positive electrode [m]

p.L_ccn = 1.2e-5;    % Thickness of negative current collector [m]
p.L_ccp = 1.6e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 5.86e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 5.22e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.75;       % Volume fraction in solid for neg. electrode 
% There's a chance that PyBaMM 0.75 eps_s = eps_s + eps_f so could split this up to be 0.58 (a value from literature) 
% Otherwise the volume fraction of filler in neg. electrode = 0
p.epsilon_s_p = 0.665;      % Volume fraction in solid for pos. electrode

% Porosity
p.epsilon_e_n = 0.25;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.47;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.335;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 1-p.epsilon_s_n-p.epsilon_e_n;  %p.epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 1-p.epsilon_s_p-p.epsilon_e_p;  %p.epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities (Guesses from typical values, not exact)
rho_sn = 1657;    % Solid phase in negative electrode [kg/m^3]
rho_sp = 3262;    % Solid phase in positive electrode [kg/m^3]
rho_ss = 397;     % Separator density [kg/m^3]
rho_ccn = 8960;   % Current collector in negative electrode
rho_ccp = 2700;   % Current collector in positive electrode

% Compute cell mass [kg/m^2]
m_n = p.L_n * rho_sn;
m_s = p.L_s * rho_ss;
m_p = p.L_p * rho_sp;
m_cc = rho_ccn*p.L_ccn + rho_ccp*p.L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;

clear rho_sn rho_sp rho_ss rho_ccn rho_ccp m_n m_s m_p m_cc;
%% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 3.3e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 4e-15;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusional conductivity in electrolyte
% f_c/a = Activity coefficient of the electrolyte
% PyBaMM: d f_c/a / c_e is either neglected or embedded in electrolyteCond.m 
% So following PyBaMM dynamics. dactivity = d ln f_ca / d ln c_e = 0
p.dactivity = 0.0;

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 215;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 0.18;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 0.2594;       % Transference number
p.Faraday = 96485.33212;    % Faraday's constant, [Coulumbs/mol]
p.Area_n = 0.0997;%0.0997;%0.1033; %0.1027; 1033        % Electrode current collector area [m^2] (Is this not electrode cross-sectional area)
p.Area_s = 0.0997;%0.0997;%0.1023; % Average of area_n and area_p
p.Area_p = 0.0997;%0.0997;%0.1033; %0.10135;        % Electrode current collector area [m^2] (Is this not electrode cross-sectional area)
p.Area_ref = mean([p.Area_n, p.Area_p]);

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 0.001;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 0;         % Contact Resistance/Current Collector Resistance, [Ohms-m^2]

% Nominal Reaction rates
p.k_n0 = 6.48e-7;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 3.42e-6; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
%% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 35e3;
p.E.kp = 17800;
p.E.Dsn = 3.03e4; % Not given by Chen et al (2020), so taken from Ecker et al. (2015) instead
p.E.Dsp = 25000; % Not given by Chen et al (2020), so taken from Ecker et al. (2015) instead
p.E.De = 17000; % Not given by Chen et al (2020), so taken from Ecker et al. (2015) instead
p.E.kappa_e = 17000; % Not given by Chen et al (2020), so taken from Ecker et al. (2015) instead

% Reference temperature
p.T_ref = 298.15; %[K]

% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7;    % [J/K]
p.C2 = 4.5;     % [J/K]
p.h12 = 10; %1.9386; % [W/K]
p.h2a = 21.45;  % [W/K]

%% Concentrations
p.c_s_n_max = 33133.0;   % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 63104.0;    % Max concentration in cathode, [mol/m^3]

p.csn0 = 29866.0;           % Initial concentration in anode, [mol/m^3]
p.csp0 = 17038.0;           % Initial concentration in cathode [mol/m^3]

p.n_Li_s = 0.2840;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
p.volt_max = 4.2;
p.volt_min = 2.5; 

%% Aging submodel parameters

%   SEI Layer Growth model
%   Adopted from Ramadass et al (2004) [Univ of South Carolina]
%   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
%   DOI: 10.1149/1.1634273
%   NOTE: These parameters have NOT been experimentally validated by eCAL

p.kappa_P = 1;      % [S/m] conductivity of side rxn product
p.M_P = 7.3e1;      % [kg/mol] molecular weight of side rxn product
p.rho_P = 2.1e3;    % [kg/m^3] mass density of side rxn product
p.i0s = 0; %1.5e-6;     % [A/m^2] exchange current density of side rxn
p.Us = 0.4;         % [V] reference potential of side rxn

%% Discretization parameters
% Discrete time step
p.delta_t = 1;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along x-coordinate
p.Nr = 30;
p.Nxn = 70;
p.Nxs = 35;
p.Nxp = 70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;


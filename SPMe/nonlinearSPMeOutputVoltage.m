%% Nonlinear output for SPMe voltage
%   Created April 28, 2025 by Wonoo Choo

function [V, varargout] = nonlinearSPMeOutputVoltage(p,c_ss_n,c_ss_p, c_ex,I)

    % Stochiometric Concentration Ratio
    theta_n = c_ss_n / p.c_s_n_max;
    theta_p = c_ss_p / p.c_s_p_max;
    
    % Equilibrium Potential
    Unref = refPotentialAnode(p,theta_n);
    Upref = refPotentialCathode(p,theta_p);
    
    % Electrolyte concentration
    % Separate and aggregate
    ce0n = c_ex(1);
    cens = c_ex(p.Nxn+1);
    cesp = c_ex(p.Nxn+p.Nxs+1);
    ce0p = c_ex(end);

    % Average electrolyte concentrations
    cen_bar = mean(c_ex(1:p.Nxn+1));
    ces_bar = mean(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1)));
    cep_bar = mean(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1)));

    % ce_bar = zeros(p.Nx,1);
    ce_bar = casadi.MX.ones(p.Nx, 1);
    ce_bar(1:p.Nxn) = cen_bar;
    ce_bar((p.Nxn+1):(p.Nxn+p.Nxs-1)) = ces_bar;
    ce_bar(p.Nxn+p.Nxs+1:end) = cep_bar;

    [i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,ce_bar,p.T_amb);
    RTaF=(p.R*p.T_amb)/(p.alph*p.Faraday);
    
    % Voltage
    eta_n = RTaF * asinh(I / (2*p.a_s_n*p.Area_n*p.L_n*i_0n(1)));
    eta_p = RTaF * asinh(-I / (2*p.a_s_p*p.Area_p*p.L_p*i_0p(end)));
    V_noVCE =  eta_p - eta_n + Upref - Unref ...
        - (p.R_f_n/(p.a_s_n*p.L_n*p.Area_n) + p.R_f_p/(p.a_s_p*p.L_p*p.Area_p))*I;
    
    % Overpotentials due to electrolyte subsystem
    kap_n = electrolyteCond(p, cen_bar, p.T_amb);
    kap_s = electrolyteCond(p, ces_bar, p.T_amb);
    kap_p = electrolyteCond(p, cep_bar, p.T_amb);
    
    kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
    kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
    kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);

    % Overpotential due to electrolyte conductivity
    V_electrolyteCond = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*I; ...
            
    % Overpotential due to electrolyte polarization
    V_electrolytePolar = (2*p.R*p.T_amb)/(p.Faraday) * ...
                            (1 * (log(cens) - log(ce0n)) ...
                            +1 * (log(cesp) - log(cens)) ...
                            +1 * (log(ce0p) - log(cesp)));

    % Add 'em up!
    V = V_noVCE + V_electrolyteCond + V_electrolytePolar;

    eta_pl = eta_n + Unref;
    if (nargout == 2)
        varargout{1} = eta_pl;
    elseif (nargout == 3)
        varargout{1} = eta_pl;
        varargout{2} = eta_n;
    elseif (nargout > 3)
        varargout{1} = eta_pl;
        varargout{2} = eta_n;
        varargout{3} = eta_p;
        varargout{4} = Unref;
        varargout{5} = Upref;
        varargout{6} = V_noVCE;
        varargout{7} = V_electrolyteCond;
        varargout{8} = V_electrolytePolar;
    end
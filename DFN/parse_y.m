function [Volt, SOC, eta_pl, c_ss_n, c_ss_p, c_avg_n, c_avg_p, c_ex, eta_n, eta_p, Unref] = parse_y(y, t_last, p)
    % Vector lengths
    Nn = p.Nxn - 1;
    Np = p.Nxp - 1;
    Nx = p.Nx - 3;

    Volt = y(1, 1:t_last);
    SOC = y(2, 1:t_last);
    eta_pl = y(3 : Nn + 2, 1:t_last);
    c_ss_n = y(Nn + 3 : 2*Nn + 2, 1:t_last);
    c_ss_p = y(2*Nn + 3 : 2*Nn + Np + 2, 1:t_last);
    
    c_avg_n = y(2*Nn + Np + 3 : 3*Nn + Np + 2, 1:t_last);
    c_avg_p = y(3*Nn + Np + 3 : 3*Nn + 2*Np + 2, 1:t_last);
    
    c_ex = y(3*Nn + 2*Np + 3 : 3*Nn + 2*Np + Nx + 6, 1:t_last);
    
    eta_n = y(3*Nn + 2*Np + Nx + 7 : 4*Nn + 2*Np + Nx + 6, 1:t_last);
    eta_p = y(4*Nn + 2*Np + Nx + 7 : 4*Nn + 3*Np + Nx + 6, 1:t_last);

    Unref = y(4*Nn + 3*Np + Nx + 7 :  5*Nn + 3*Np + Nx + 6, 1:t_last);
end
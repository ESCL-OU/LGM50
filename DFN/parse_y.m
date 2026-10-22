function [Volt, SOC, eta_pl, c_ss_n, c_ss_p, c_avg_n, c_avg_p, c_ex, eta_n, eta_p, Unref, varargout] = parse_y(y, t_last, p)
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

    if nargout > 11
        idx = 5*Nn + 3*Np + Nx + 7;
        deg = struct();
        if size(y, 1) >= idx + 10*Nn + 3
            deg.LDP_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.j_int_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.j_sei_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.j_lpl_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.j_tot_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.delta_sei_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.delta_li_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.epsilon_e_n_dyn = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.phi_e_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.phi_s_n = y(idx:idx+Nn-1, 1:t_last); idx = idx + Nn;
            deg.LDP_sep_n = y(idx, 1:t_last); idx = idx + 1;
            deg.j_lpl_sep_n = y(idx, 1:t_last); idx = idx + 1;
            deg.delta_li_sep_n = y(idx, 1:t_last); idx = idx + 1;
            deg.epsilon_e_sep_n = y(idx, 1:t_last);
        else
            deg.LDP_n = zeros(Nn, t_last);
            deg.j_int_n = zeros(Nn, t_last);
            deg.j_sei_n = zeros(Nn, t_last);
            deg.j_lpl_n = zeros(Nn, t_last);
            deg.j_tot_n = zeros(Nn, t_last);
            deg.delta_sei_n = zeros(Nn, t_last);
            deg.delta_li_n = zeros(Nn, t_last);
            deg.epsilon_e_n_dyn = p.epsilon_e_n * ones(Nn, t_last);
            deg.phi_e_n = zeros(Nn, t_last);
            deg.phi_s_n = zeros(Nn, t_last);
            deg.LDP_sep_n = zeros(1, t_last);
            deg.j_lpl_sep_n = zeros(1, t_last);
            deg.delta_li_sep_n = zeros(1, t_last);
            deg.epsilon_e_sep_n = p.epsilon_e_n * ones(1, t_last);
        end
        varargout{1} = deg;
    end
end
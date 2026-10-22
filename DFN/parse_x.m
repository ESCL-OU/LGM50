function [c_s_n, c_s_p, c_e, T, varargout] = parse_x(x, t_last, p)
    % Vector lengths
    Ncsn = p.PadeOrder * (p.Nxn-1);
    Ncsp = p.PadeOrder * (p.Nxp-1);
    Nce = p.Nx - 3;
    Nc = Ncsn+Ncsp+Nce;
    Nn = p.Nxn - 1;
    use_degradation = isfield(p, 'dfn_degradation_enabled') && p.dfn_degradation_enabled;
    % Parse out States
    c_s_n = x(1:Ncsn, 1:t_last);
    c_s_p = x(Ncsn+1:Ncsn+Ncsp, 1:t_last);
    c_e = x(Ncsn+Ncsp+1:Nc, 1:t_last);
    if use_degradation
        deg.delta_sei_n = x(Nc+1:Nc+Nn, 1:t_last);
        deg.delta_li_n = x(Nc+Nn+1:Nc+2*Nn, 1:t_last);
        deg.epsilon_e_n_dyn = x(Nc+2*Nn+1:Nc+3*Nn, 1:t_last);
        T = x(Nc+3*Nn+1, 1:t_last);
    else
        deg.delta_sei_n = zeros(Nn, t_last);
        deg.delta_li_n = zeros(Nn, t_last);
        deg.epsilon_e_n_dyn = p.epsilon_e_n * ones(Nn, t_last);
        T = x(end, 1:t_last);
    end
    if nargout > 4
        varargout{1} = deg;
    end
end
function [c_s_n, c_s_p, c_e, T] = parse_x(x, t_last, p)
    % Vector lengths
    Ncsn = p.PadeOrder * (p.Nxn-1);
    Ncsp = p.PadeOrder * (p.Nxp-1);
    Nce = p.Nx - 3;
    Nc = Ncsn+Ncsp+Nce;
    % Parse out States
    c_s_n = x(1:Ncsn, 1:t_last);
    c_s_p = x(Ncsn+1:Ncsn+Ncsp, 1:t_last);
    c_e = x(Ncsn+Ncsp+1:Nc, 1:t_last);
    T = x(end, 1:t_last);
end
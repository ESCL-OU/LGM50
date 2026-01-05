function [phi_s_n, phi_s_p, i_en, i_ep, phi_e, jn, jp] = parse_z(z, t_last, p)
    % Vector lengths
    Nn = p.Nxn - 1;
    Np = p.Nxp - 1;
    Nnp = Nn+Np;
    Nx = p.Nx - 3;

    phi_s_n = z(1:Nn, 1:t_last);
    phi_s_p = z(Nn+1:Nnp, 1:t_last);
    i_en = z(Nnp+1:Nnp+Nn, 1:t_last);
    i_ep = z(Nnp+Nn+1:2*Nnp, 1:t_last);
    phi_e = z(2*Nnp+1:2*Nnp+Nx+2, 1:t_last);
    jn = z(2*Nnp+Nx+3:2*Nnp+Nx+Nn+2, 1:t_last);
    jp = z(2*Nnp+Nx+Nn+3:end, 1:t_last);
end

function [c_e_dot, varargout] = ode_electrolyte(c_e,jn,jp,p)

    [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

    % Compute Boundary Conditions
    c_e_bcs = C_ce * c_e;

    % Separate and aggregate
    c_en = c_e(1:(p.Nxn-1));
    c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
    c_ep = c_e((p.Nxn-1)+p.Nxs : end);
    c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

    % Compute Electrolyte Diffusion Coefficient and Derivative
    [D_en,dD_en] = electrolyteDe(p, c_en, p.T_amb);
    [D_es,dD_es] = electrolyteDe(p, c_es, p.T_amb);
    [D_ep,dD_ep] = electrolyteDe(p, c_ep, p.T_amb);
    
    % ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
    D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
    dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

    D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
    dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

    D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
    dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);
    
    % Compute derivative
    c_en_dot = dD_en_eff.*(M1n*c_en + M2n*c_e_bcs(1:2)).^2 ...
        + D_en_eff.*(M3n*c_en + M4n*c_e_bcs(1:2)) + diag(M5n)*jn;

    c_es_dot = dD_es_eff.*(M1s*c_es + M2s*c_e_bcs(2:3)).^2 ...
        + D_es_eff.*(M3s*c_es + M4s*c_e_bcs(2:3));

    c_ep_dot = dD_ep_eff.*(M1p*c_ep + M2p*c_e_bcs(3:4)).^2 ...
        + D_ep_eff.*(M3p*c_ep + M4p*c_e_bcs(3:4)) + diag(M5p)*jp;

    % Assemble c_e_dot
    c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];
    
    if nargout == 2
        varargout{1} = c_ex;
    end
end
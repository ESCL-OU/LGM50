function [eta_n, eta_p, Unref, Upref, IR_n, IR_p, phi_e_n, phi_e_p, V_noVCE, V_electrolyteCond, V_electrolytePolar] = potentials(x, z, I, p)
    % Output data
    [~, ~, ~, eta_n, eta_p, Unref, Upref, IR_n, IR_p, phi_e_n, phi_e_p, V_noVCE, V_electrolyteCond, V_electrolytePolar] = dae_dfn(x, z, I, p);

end
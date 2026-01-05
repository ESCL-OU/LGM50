function [p, ecm_p, ocv_p] = load_ECM_params(z0, n_order)
    %% Electrochemical Model Parameters for Q
    run param\params_LGM50
    %% OCV curve Polynomial coefficients
    load("param\OCVdata.mat", "ocv_p");
    %% Equivalent Circuit Model Parameters;
    load("param\ECM_Params_order" + string(n_order) + ".mat", "ecm_p")

    %% Compute initial States
    z = min(max(z0, 0), 1);
    p.x0 = zeros(n_order+1, 1);
    p.x0(1) = z;
    p.y0 = polyval(ocv_p, z);
end
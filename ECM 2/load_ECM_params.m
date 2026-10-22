function [p, ecm_p, ocv_p] = load_ECM_params(z0, n_order)
    %% Electrochemical Model Parameters for Q
    base_dir = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(base_dir, 'param', 'params_LGM50.m'))
    %% OCV curve Polynomial coefficients
    load(fullfile(base_dir, 'param', 'OCVdata.mat'), "ocv_p");
    %% Equivalent Circuit Model Parameters;
    load(fullfile(base_dir, 'param', "ECM_Params_order" + string(n_order) + ".mat"), "ecm_p")

    %% Compute initial States
    z = min(max(z0, 0), 1);
    p.x0 = zeros(n_order+1, 1);
    p.x0(1) = z;
    p.y0 = polyval(ocv_p, z);
end
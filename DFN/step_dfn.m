function [x, z, y, varargout] = step_dfn(x0, z0, Cur_vec, p, verbosity)
    % Simulate DFN plant
    % Solid concentration matrices
    T = x0(end);
    [p.A_csn, p.B_csn, p.A_csp, p.B_csp, p.C_csn, p.C_csp] = c_s_mats(p, T);

    % Step-forward in time
    [x, z, stats] = cn_dfn(x0, z0, Cur_vec, p, verbosity);

    % Output data
    [~, ~, y] = dae_dfn(x, z, Cur_vec(2), p);

    if nargout == 4
        varargout{1} = stats;
    end
end
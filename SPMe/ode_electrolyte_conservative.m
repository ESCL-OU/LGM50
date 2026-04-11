function [c_e_dot, aux] = ode_electrolyte_conservative(c_e, jn, jp, p)
%ODE_ELECTROLYTE_CONSERVATIVE
% Conservative electrolyte update for the SPMe.
%
% State layout is kept the same as the current code:
%   c_e = [c_en ; c_es ; c_ep]
% where
%   c_en : (p.Nxn-1) anode interior electrolyte concentrations
%   c_es : (p.Nxs-1) separator interior electrolyte concentrations
%   c_ep : (p.Nxp-1) cathode interior electrolyte concentrations
%
% Governing semi-discrete form:
%   eps_i * dx_i * dc_i/dt = F_{i-1/2} - F_{i+1/2} + S_i * dx_i
%
% with one shared face flux per face:
%   F_{i+1/2} = - D_eff_face * (c_{i+1} - c_i) / dxc_face
%
% This guarantees discrete conservation of electrolyte lithium at the
% semi-discrete level, up to roundoff, provided the source terms are
% globally balanced and the outer boundary fluxes are zero.
%
% Notes
% - This only rewrites the electrolyte operator.
% - Solid-particle dynamics can stay exactly as they are.
% - Uses the same D_eff convention as the current code:
%       D_eff = D_e * epsilon_e^(brug - 1)
%   matching the current implementation. :contentReference[oaicite:3]{index=3}
%
% Optional second output "aux" returns diagnostic data.

    import casadi.*

    c_e = c_e(:);

    % ----------------------------
    % Sizes and region extraction
    % ----------------------------
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;

    idx_n = 1:Nn;
    idx_s = (Nn+1):(Nn+Ns);
    idx_p = (Nn+Ns+1):(Nn+Ns+Np);

    c_en = c_e(idx_n);
    c_es = c_e(idx_s);
    c_ep = c_e(idx_p);

    % ----------------------------
    % Region spacings
    % ----------------------------
    dx_n = p.L_n * p.delta_x_n;
    dx_s = p.L_s * p.delta_x_s;
    dx_p = p.L_p * p.delta_x_p;

    % ----------------------------
    % Effective diffusivities
    % ----------------------------
    [D_en, ~] = electrolyteDe(p, c_en, p.T_amb);
    [D_es, ~] = electrolyteDe(p, c_es, p.T_amb);
    [D_ep, ~] = electrolyteDe(p, c_ep, p.T_amb);

    D_en_eff = D_en .* p.epsilon_e_n.^(p.brug - 1);
    D_es_eff = D_es .* p.epsilon_e_s.^(p.brug - 1);
    D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug - 1);

    % ----------------------------
    % Source terms  S/epsilon
    % PDE form used in the current code is:
    %   dc/dt = ... + (1-t_plus)*a_s/epsilon_e * j
    % matching M5n/M5p in c_e_mats.m. :contentReference[oaicite:4]{index=4}
    % ----------------------------
    src_n = ((1 - p.t_plus) * p.a_s_n / p.epsilon_e_n) * jn * ones_like(c_en);
    src_s = zeros_like(c_es);
    src_p = ((1 - p.t_plus) * p.a_s_p / p.epsilon_e_p) * jp * ones_like(c_ep);

    % ----------------------------
    % Face fluxes
    % Sign convention:
    %   F > 0 means flux in +x direction
    %   Update uses: eps*dx*dc/dt = F_left - F_right + S*dx
    %
    % Outer boundaries are zero-flux:
    %   F_left_boundary  = 0
    %   F_right_boundary = 0
    % ----------------------------

    % Faces inside anode region
    F_n = zeros_like_len(max(Nn-1,0), c_e);
    for k = 1:(Nn-1)
        Df = harmonic_mean(D_en_eff(k), D_en_eff(k+1));
        dxc = dx_n;
        F_n(k) = -Df * (c_en(k+1) - c_en(k)) / dxc;
    end

    % Faces inside separator region
    F_s = zeros_like_len(max(Ns-1,0), c_e);
    for k = 1:(Ns-1)
        Df = harmonic_mean(D_es_eff(k), D_es_eff(k+1));
        dxc = dx_s;
        F_s(k) = -Df * (c_es(k+1) - c_es(k)) / dxc;
    end

    % Faces inside cathode region
    F_p = zeros_like_len(max(Np-1,0), c_e);
    for k = 1:(Np-1)
        Df = harmonic_mean(D_ep_eff(k), D_ep_eff(k+1));
        dxc = dx_p;
        F_p(k) = -Df * (c_ep(k+1) - c_ep(k)) / dxc;
    end

    % Shared interface face: anode | separator
    Df_ns = harmonic_mean(D_en_eff(end), D_es_eff(1));
    dxc_ns = 0.5 * dx_n + 0.5 * dx_s;
    F_ns = -Df_ns * (c_es(1) - c_en(end)) / dxc_ns;

    % Shared interface face: separator | cathode
    Df_sp = harmonic_mean(D_es_eff(end), D_ep_eff(1));
    dxc_sp = 0.5 * dx_s + 0.5 * dx_p;
    F_sp = -Df_sp * (c_ep(1) - c_es(end)) / dxc_sp;

    % Outer boundary fluxes
    F_left_bc  = scalar_zero_like(c_e);
    F_right_bc = scalar_zero_like(c_e);

    % ----------------------------
    % Conservative updates
    % eps*dx*dc/dt = F_left - F_right + S*dx
    % ----------------------------
    c_en_dot = zeros_like(c_en);
    c_es_dot = zeros_like(c_es);
    c_ep_dot = zeros_like(c_ep);

    % Anode cells
    for i = 1:Nn
        if i == 1
            F_left = F_left_bc;
        else
            F_left = F_n(i-1);
        end

        if i == Nn
            F_right = F_ns;
        else
            F_right = F_n(i);
        end

        c_en_dot(i) = (F_left - F_right) / (p.epsilon_e_n * dx_n) + src_n(i);
    end

    % Separator cells
    for i = 1:Ns
        if i == 1
            F_left = F_ns;
        else
            F_left = F_s(i-1);
        end

        if i == Ns
            F_right = F_sp;
        else
            F_right = F_s(i);
        end

        c_es_dot(i) = (F_left - F_right) / (p.epsilon_e_s * dx_s) + src_s(i);
    end

    % Cathode cells
    for i = 1:Np
        if i == 1
            F_left = F_sp;
        else
            F_left = F_p(i-1);
        end

        if i == Np
            F_right = F_right_bc;
        else
            F_right = F_p(i);
        end

        c_ep_dot(i) = (F_left - F_right) / (p.epsilon_e_p * dx_p) + src_p(i);
    end

    c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

    % Optional diagnostics
    if nargout > 1
        aux = struct();
        aux.c_en = c_en;
        aux.c_es = c_es;
        aux.c_ep = c_ep;
        aux.D_en_eff = D_en_eff;
        aux.D_es_eff = D_es_eff;
        aux.D_ep_eff = D_ep_eff;
        aux.F_n = F_n;
        aux.F_s = F_s;
        aux.F_p = F_p;
        aux.F_ns = F_ns;
        aux.F_sp = F_sp;
        aux.src_n = src_n;
        aux.src_s = src_s;
        aux.src_p = src_p;
    end
end

% ============================================================
% Helpers
% ============================================================

function y = harmonic_mean(a, b)
% Harmonic mean, suitable for a shared diffusive face.
% For positive diffusivities, this is robust across material jumps.
    y = 2 .* a .* b ./ (a + b);
end

function z = zeros_like(x)
    import casadi.*
    if isa(x, 'casadi.MX')
        z = MX.zeros(size(x,1), size(x,2));
    elseif isa(x, 'casadi.SX')
        z = SX.zeros(size(x,1), size(x,2));
    else
        z = zeros(size(x));
    end
end

function z = zeros_like_len(n, ref)
    import casadi.*
    if isa(ref, 'casadi.MX')
        z = MX.zeros(n, 1);
    elseif isa(ref, 'casadi.SX')
        z = SX.zeros(n, 1);
    else
        z = zeros(n, 1);
    end
end

function z = scalar_zero_like(ref)
    import casadi.*
    if isa(ref, 'casadi.MX')
        z = MX.zeros(1,1);
    elseif isa(ref, 'casadi.SX')
        z = SX.zeros(1,1);
    else
        z = 0;
    end
end

function y = ones_like(x)
    import casadi.*
    if isa(x, 'casadi.MX')
        y = MX.ones(size(x,1), size(x,2));
    elseif isa(x, 'casadi.SX')
        y = SX.ones(size(x,1), size(x,2));
    else
        y = ones(size(x));
    end
end
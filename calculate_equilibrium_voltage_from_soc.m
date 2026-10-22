function out = calculate_equilibrium_voltage_from_soc(target_soc)
%CALCULATE_EQUILIBRIUM_VOLTAGE_FROM_SOC Single-SoC equilibrium voltage.
%   out = calculate_equilibrium_voltage_from_soc(target_soc)
%
%   Procedure:
%   1) Find endpoint concentrations from init_cs at Vmin (SoC0) and
%      Vmax (SoC100).
%   2) Linearly interpolate anode/cathode concentrations at target SoC.
%   3) Evaluate equilibrium voltage using electrode reference potentials.
%
%   Output fields:
%   - SOC_target, Veq
%   - cn_soc0, cp_soc0, cn_soc100, cp_soc100
%   - cn_interp, cp_interp, theta_n, theta_p, U_anode, U_cathode
%   - p (parameter struct used for calculation)

    if nargin < 1
        error('target_soc input is required.');
    end
    if ~isscalar(target_soc)
        error('target_soc must be a single scalar value.');
    end
    if target_soc < 0 || target_soc > 1
        error('target_soc must be in [0, 1].');
    end

    base_dir = fileparts(mfilename('fullpath'));
    param_dir = fullfile(base_dir, 'param');
    if ~isfolder(param_dir)
        error('param folder not found at %s', param_dir);
    end
    addpath(param_dir);

    run(fullfile(param_dir, 'params_LGM50.m'));

    % Endpoint concentrations from voltage limits.
    [cn_soc0, cp_soc0] = init_cs(p, p.volt_min);
    [cn_soc100, cp_soc100] = init_cs(p, p.volt_max);

    % Same interpolation procedure for anode and cathode.
    cn_interp = cn_soc0 + target_soc * (cn_soc100 - cn_soc0);
    cp_interp = cp_soc0 + target_soc * (cp_soc100 - cp_soc0);

    theta_n = cn_interp / p.c_s_n_max;
    theta_p = cp_interp / p.c_s_p_max;

    U_anode = refPotentialAnode(p, theta_n);
    U_cathode = refPotentialCathode(p, theta_p);
    Veq = U_cathode - U_anode;

    out = struct();
    out.SOC_target = target_soc;
    out.Veq = Veq;
    out.cn_soc0 = cn_soc0;
    out.cp_soc0 = cp_soc0;
    out.cn_soc100 = cn_soc100;
    out.cp_soc100 = cp_soc100;
    out.cn_interp = cn_interp;
    out.cp_interp = cp_interp;
    out.theta_n = theta_n;
    out.theta_p = theta_p;
    out.U_anode = U_anode;
    out.U_cathode = U_cathode;
    out.p = p;
end

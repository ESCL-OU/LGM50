function out = test_dfn_electrolyte_li_inventory_modes()
% Diagnose electrolyte Li-ion inventory drift / negative c_e in DFN.
%
% Compares:
%   1. no SEI + no Li plating
%   2. SEI only
%   3. Li plating only
%   4. SEI + Li plating
%
% Uses DFN + padded_analytic Jacobian.

init_soc = 0.1;
profile = repmat({'CC', -1.0, 'CC', 1.0}, 1, 5);

delta_t = 1;
t_amb = 25;

Nr = 10;
Nxn = 10;
Nxs = 5;
Nxp = 10;

modes = { ...
    'no_sei_no_lpl', false, false; ...
    'sei_only',      false, false; ... % handled by safari + no plating
    'lpl_only',      true,  true;  ...
    'sei_lpl',       true,  true};

out = struct();

for m = 1:size(modes,1)
    name = modes{m,1};

    switch name
        case 'no_sei_no_lpl'
            degradation_model = 'none';
            enable_lpl = false;
            enable_porosity = false;

        case 'sei_only'
            degradation_model = 'safari';
            enable_lpl = false;
            enable_porosity = true;

        case 'lpl_only'
            % NOTE: this assumes your code can disable SEI separately.
            % If not, this mode will still include Safari SEI.
            degradation_model = 'safari';
            enable_lpl = true;
            enable_porosity = true;

        case 'sei_lpl'
            degradation_model = 'safari';
            enable_lpl = true;
            enable_porosity = true;
    end

    fprintf('\nRunning mode: %s\n', name);
    tic;

    sim = sim_profile_sequence_LGM50('DFN', init_soc, profile, ...
        'Repeat', 1, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', degradation_model, ...
        'EnableLithiumPlating', enable_lpl, ...
        'EnableAnodePorosityLoss', enable_porosity, ...
        'DFNLplCouplingMode', 'spme_compatible', ...
        'DFNJacobianMode', 'padded_analytic', ...
        'MemoryLight', false, ...
        'ShowProgress', true, ...
        'SaveOutput', false);

    elapsed = toc;
    diag = compute_electrolyte_inventory_diagnostics(sim);

    out.(name).sim = sim;
    out.(name).diag = diag;
    out.(name).elapsed_s = elapsed;

    fprintf('Finished %s in %.2f s | min c_e = %.3e | electrolyte Li drift = %.3e mol\n', ...
        name, elapsed, min(diag.c_e_min), diag.n_electrolyte_mol(end)-diag.n_electrolyte_mol(1));
end

plot_electrolyte_inventory_diagnostics(out);

end

function diag = compute_electrolyte_inventory_diagnostics(sim)

    p = sim.p;
    c_e = sim.c_e;
    t = sim.time(:);
    
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    
    Ln = p.L_n;
    Ls = p.L_s;
    Lp = p.L_p;
    if isfield(p, 'Area')
        A = p.Area;
    elseif isfield(p, 'Area_ref')
        A = p.Area_ref;
    else
        warning('No p.Area or p.Area_ref found. Using A = 1, so inventory is per unit area.');
        A = 1;
    end
    
    dx_n = Ln / Nn;
    dx_s = Ls / Ns;
    dx_p = Lp / Np;
    
    vol_weights = [ ...
        p.epsilon_e_n * A * dx_n * ones(Nn,1); ...
        p.epsilon_e_s * A * dx_s * ones(Ns,1); ...
        p.epsilon_e_p * A * dx_p * ones(Np,1)];
    
    if isfield(sim, 'dfn_outputs') && isfield(sim.dfn_outputs, 'epsilon_e_n_dyn')
        
        eps_n = [];

        if isfield(sim, 'x') && ~isempty(sim.x)
            try
                [~, ~, ~, ~, xdeg] = parse_x(sim.x, numel(t), p);
                if isfield(xdeg, 'epsilon_e_n_dyn')
                    eps_n = xdeg.epsilon_e_n_dyn;
                end
            catch
                eps_n = [];
            end
        end
        
        if isempty(eps_n) && isfield(sim, 'dfn_outputs') && isfield(sim.dfn_outputs, 'epsilon_e_n_dyn')
            tmp = sim.dfn_outputs.epsilon_e_n_dyn;
            if size(tmp,2) == numel(t)
                eps_n = tmp;
            end
        end
        
        if ~isempty(eps_n) && size(eps_n,2) == numel(t)
            vol_weights_time = zeros(size(c_e));
            for k = 1:numel(t)
                vol_weights_time(:,k) = [ ...
                    eps_n(:,k) * A * dx_n; ...
                    p.epsilon_e_s * A * dx_s * ones(Ns,1); ...
                    p.epsilon_e_p * A * dx_p * ones(Np,1)];
            end
        else
            vol_weights_time = repmat(vol_weights, 1, numel(t));
        end
        
        vol_weights_time = zeros(size(c_e));
        for k = 1:numel(t)
            vol_weights_time(:,k) = [ ...
                eps_n(:,k) * A * dx_n; ...
                p.epsilon_e_s * A * dx_s * ones(Ns,1); ...
                p.epsilon_e_p * A * dx_p * ones(Np,1)];
        end
    else
        vol_weights_time = repmat(vol_weights, 1, numel(t));
    end
    
    n_electrolyte_mol = sum(c_e .* vol_weights_time, 1).';
    
    diag = struct();
    diag.time_s = t;
    diag.time_h = t/3600;
    diag.n_electrolyte_mol = n_electrolyte_mol;
    diag.n_electrolyte_mol_drift = n_electrolyte_mol - n_electrolyte_mol(1);
    diag.c_e_min = min(c_e, [], 1).';
    diag.c_e_max = max(c_e, [], 1).';
    diag.c_e_mean = mean(c_e, 1).';
    diag.voltage = sim.voltage(:);
    diag.current = sim.current(:);
    diag.soc = sim.soc(:);
    
    if isfield(sim, 'dfn_outputs')
        d = sim.dfn_outputs;
        maybe = {'LDP_sep_n','j_lpl_sep_n','delta_li_sep_n','epsilon_e_sep_n'};
        for i = 1:numel(maybe)
            fn = maybe{i};
            if isfield(d, fn)
                diag.(fn) = d.(fn)(:);
            end
        end
    end

end

function plot_electrolyte_inventory_diagnostics(out)

    names = fieldnames(out);
    
    figure('Name','DFN electrolyte Li-ion inventory diagnostics','Color','w');
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
    
    for i = 1:numel(names)
        name = names{i};
        d = out.(name).diag;
    
        nexttile(1); hold on; grid on; box on;
        plot(d.time_h, d.n_electrolyte_mol_drift, 'DisplayName', name);
        ylabel('\Delta n_{e,Li+} [mol]');
        title('Electrolyte Li inventory drift');
    
        nexttile(2); hold on; grid on; box on;
        plot(d.time_h, d.c_e_min, 'DisplayName', name);
        yline(0,'k--','HandleVisibility','off');
        ylabel('min c_e [mol/m^3]');
        title('Minimum electrolyte concentration');
    
        nexttile(3); hold on; grid on; box on;
        plot(d.time_h, d.c_e_mean, 'DisplayName', name);
        ylabel('mean c_e [mol/m^3]');
        title('Mean electrolyte concentration');
    
        nexttile(4); hold on; grid on; box on;
        plot(d.time_h, d.voltage, 'DisplayName', name);
        ylabel('Voltage [V]');
        title('Voltage');
    
        nexttile(5); hold on; grid on; box on;
        plot(d.time_h, d.current, 'DisplayName', name);
        ylabel('Current [A]');
        xlabel('Time [h]');
        title('Current');
    
        nexttile(6); hold on; grid on; box on;
        plot(d.time_h, d.soc, 'DisplayName', name);
        ylabel('SOC [-]');
        xlabel('Time [h]');
        title('SOC');
    end
    
    for k = 1:6
        nexttile(k);
        legend('Location','best','Interpreter','none');
    end
    
    figure('Name','DFN degradation diagnostics vs electrolyte inventory','Color','w');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    
    for i = 1:numel(names)
        name = names{i};
        d = out.(name).diag;
    
        if isfield(d,'LDP_sep_n')
            nexttile(1); hold on; grid on; box on;
            plot(d.time_h, d.LDP_sep_n, 'DisplayName', name);
            ylabel('LDP_{sep} [V]');
            title('Separator-side plating potential');
        end
    
        if isfield(d,'j_lpl_sep_n')
            nexttile(2); hold on; grid on; box on;
            semilogy(d.time_h, max(abs(d.j_lpl_sep_n), realmin), 'DisplayName', name);
            ylabel('|j_{lpl,sep}|');
            title('Separator-side plating flux');
        end
    
        if isfield(d,'delta_li_sep_n')
            nexttile(3); hold on; grid on; box on;
            plot(d.time_h, d.delta_li_sep_n*1e9, 'DisplayName', name);
            ylabel('\delta_{Li,sep} [nm]');
            xlabel('Time [h]');
            title('Plated Li thickness');
        end
    
        if isfield(d,'epsilon_e_sep_n')
            nexttile(4); hold on; grid on; box on;
            plot(d.time_h, d.epsilon_e_sep_n, 'DisplayName', name);
            ylabel('\epsilon_{e,n,sep}');
            xlabel('Time [h]');
            title('Separator-side anode porosity');
        end
    end
    
    for k = 1:4
        nexttile(k);
        legend('Location','best','Interpreter','none');
    end

end
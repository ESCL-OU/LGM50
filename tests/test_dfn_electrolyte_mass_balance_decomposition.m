function out = test_dfn_electrolyte_mass_balance_decomposition()
% Diagnose why DFN electrolyte concentration mean/inventory drifts.

init_soc = 0.1;
profile = {'CC', -1.0, 'CC', 1.0};

sim = sim_profile_sequence_LGM50('DFN', init_soc, profile, ...
    'Repeat', 1, ...
    'DeltaT', 1, ...
    'Tamb', 25, ...
    'Nr', 10, 'Nxn', 10, 'Nxs', 5, 'Nxp', 10, ...
    'DegradationModel', 'none', ...
    'EnableLithiumPlating', false, ...
    'EnableAnodePorosityLoss', false, ...
    'DFNJacobianMode', 'padded_analytic', ...
    'MemoryLight', false, ...
    'ShowProgress', true, ...
    'SaveOutput', false);

out = decompose_electrolyte_balance(sim);
plot_electrolyte_balance_decomposition(out);

end

function out = decompose_electrolyte_balance(sim)

p = sim.p;
t = sim.time(:);
c_e = sim.c_e;
z = sim.z;

Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nce = p.Nx - 3;

if isfield(p,'Area_ref')
    A = p.Area_ref;
else
    A = 1;
end

dx_n = p.L_n / Nn;
dx_s = p.L_s / Ns;
dx_p = p.L_p / Np;

w = [ ...
    p.epsilon_e_n * A * dx_n * ones(Nn,1); ...
    p.epsilon_e_s * A * dx_s * ones(Ns,1); ...
    p.epsilon_e_p * A * dx_p * ones(Np,1)];

n_e = sum(c_e .* w, 1).';
mean_ce = mean(c_e,1).';
sum_ce = sum(c_e,1).';
min_ce = min(c_e,[],1).';
max_ce = max(c_e,[],1).';

int_f_total = nan(numel(t),1);
int_f_diff = nan(numel(t),1);
int_f_source_n = nan(numel(t),1);
int_f_source_p = nan(numel(t),1);
int_f_source_total = nan(numel(t),1);

int_source_raw_n = nan(numel(t),1);
int_source_raw_p = nan(numel(t),1);
int_source_raw_total = nan(numel(t),1);

for k = 1:numel(t)
    xk = sim.x(:,k);
    zk = sim.z(:,k);
    Cur = sim.current(k);

    [f,~] = dae_dfn(xk, zk, Cur, p);

    Ncsn = p.PadeOrder * (p.Nxn-1);
    Ncsp = p.PadeOrder * (p.Nxp-1);
    ind_ce = Ncsn + Ncsp + (1:Nce);

    ce_dot_total = f(ind_ce);
    int_f_total(k) = sum(w .* ce_dot_total);

    [~,~,~,~,~,jn,jp] = parse_z(zk,1,p);

    source_n = p.ce.M5n * jn;
    source_p = p.ce.M5p * jp;
    
    source_full = [source_n; zeros(Ns,1); source_p];
    
    % Physical mol/s source after multiplying by epsilon*volume.
    % Since M5 already contains 1/epsilon_e, use w .* source_full.
    int_f_source_n(k) = sum(w(1:Nn) .* source_n);
    int_f_source_p(k) = sum(w(Nn+Ns+1:end) .* source_p);
    int_f_source_total(k) = int_f_source_n(k) + int_f_source_p(k);
    
    % Diffusion/boundary contribution to inventory.
    ce_dot_diff = ce_dot_total - source_full;
    int_f_diff(k) = sum(w .* ce_dot_diff);
    
    % Additional unweighted diagnostic: should NOT be interpreted as mol/s,
    % but helps identify whether source rows themselves are imbalanced.
    int_source_raw_n(k) = sum(source_n);
    int_source_raw_p(k) = sum(source_p);
    int_source_raw_total(k) = int_source_raw_n(k) + int_source_raw_p(k);
end

dn_dt_numeric = [nan; diff(n_e)./diff(t)];
accum_total = cumtrapz(t, int_f_total);
accum_diff  = cumtrapz(t, int_f_diff);
accum_src   = cumtrapz(t, int_f_source_total);
accum_src_n = cumtrapz(t, int_f_source_n);
accum_src_p = cumtrapz(t, int_f_source_p);

out = struct();
out.sim = sim;
out.time_s = t;
out.time_h = t/3600;
out.c_e = c_e;
out.n_e = n_e;
out.n_e_drift = n_e;
out.mean_ce = mean_ce;
out.sum_ce = sum_ce;
out.min_ce = min_ce;
out.max_ce = max_ce;
out.dn_dt_numeric = dn_dt_numeric;
out.int_f_total = int_f_total;
out.int_f_diff = int_f_diff;
out.int_f_source_n = int_f_source_n;
out.int_f_source_p = int_f_source_p;
out.int_f_source_total = int_f_source_total;
out.accum_total = accum_total;
out.accum_diff = accum_diff;
out.accum_src = accum_src;
out.accum_src_n = accum_src_n;
out.accum_src_p = accum_src_p;
out.int_source_raw_n = int_source_raw_n;
out.int_source_raw_p = int_source_raw_p;
out.int_source_raw_total = int_source_raw_total;

fprintf('\nElectrolyte diagnostic summary:\n');
fprintf('Initial mean c_e: %.6g\n', mean_ce(1));
fprintf('Final mean c_e:   %.6g\n', mean_ce(end));
fprintf('Initial n_e:      %.6e mol\n', n_e(1));
fprintf('Final n_e:        %.6e mol\n', n_e(end));
fprintf('Drift n_e:        %.6e mol\n', n_e(end)-n_e(1));
fprintf('Max |int diffusion contribution|: %.6e mol/s\n', max(abs(int_f_diff)));
fprintf('Max |int source n|:              %.6e mol/s\n', max(abs(int_f_source_n)));
fprintf('Max |int source p|:              %.6e mol/s\n', max(abs(int_f_source_p)));
fprintf('Max |int source total|:          %.6e mol/s\n', max(abs(int_f_source_total)));

end

function plot_electrolyte_balance_decomposition(out)

t = out.time_h;

figure('Name','Electrolyte mass balance decomposition','Color','w');
tiledlayout(4,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on; box on;
plot(t, out.mean_ce);
ylabel('mean c_e [mol/m^3]');
title('Mean electrolyte concentration');

nexttile; hold on; grid on; box on;
plot(t, out.n_e_drift);
ylabel('\Delta n_e [mol]');
title('Electrolyte Li-ion inventory drift');

nexttile; hold on; grid on; box on;
plot(t, out.min_ce, 'DisplayName','min');
plot(t, out.max_ce, 'DisplayName','max');
ylabel('c_e [mol/m^3]');
title('Min/max c_e');
legend;

nexttile; hold on; grid on; box on;
plot(t, out.dn_dt_numeric, 'DisplayName','numeric dN/dt');
plot(t, out.int_f_total, '--', 'DisplayName','integral f_{ce}');
ylabel('mol/s');
title('Inventory derivative check');
legend;

nexttile; hold on; grid on; box on;
plot(t, out.int_f_diff, 'DisplayName','diffusion part');
plot(t, out.int_f_source_total, 'DisplayName','source total');
plot(t, out.int_f_source_n, '--', 'DisplayName','source anode');
plot(t, out.int_f_source_p, '--', 'DisplayName','source cathode');
ylabel('mol/s');
xlabel('Time [h]');
title('Contribution decomposition');
legend;

nexttile; hold on; box on;
imagesc(t, 1:size(out.c_e,1), out.c_e);
set(gca,'YDir','normal');
xlabel('Time [h]');
ylabel('c_e node');
title('c_e heatmap');
colorbar;

nexttile; hold on; grid on; box on;
plot(t, out.n_e - out.n_e(1), 'k-', 'LineWidth', 1.5, ...
    'DisplayName','actual inventory drift');
plot(t, out.accum_total - out.accum_total(1), '--', 'LineWidth', 1.2, ...
    'DisplayName','accum total f_{ce}');
plot(t, out.accum_diff - out.accum_diff(1), 'LineWidth', 1.2, ...
    'DisplayName','accum diffusion/boundary');
plot(t, out.accum_src - out.accum_src(1), 'LineWidth', 1.2, ...
    'DisplayName','accum source total');
ylabel('Accumulated drift [mol]');
xlabel('Time [h]');
title('Accumulated source of electrolyte inventory drift');
legend('Location','best','Interpreter','none');

end
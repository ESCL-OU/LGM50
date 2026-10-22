function out = test_dfn_jacobian_mode_speed_accuracy()
% Compare DFN degradation Jacobian modes:
%   finite_difference vs padded_analytic
%
% Purpose:
%   Check whether padded_analytic gives similar trajectories with much faster runtime.

init_soc = 0.1;
profile = {'CC', -1.0};   % short charge segment first
delta_t = 1;
t_amb = 25;

Nr = 8;
Nxn = 6;
Nxs = 3;
Nxp = 6;

modes = {'finite_difference', 'padded_analytic'};
out = struct();

for i = 1:numel(modes)
    mode = modes{i};

    fprintf('\nRunning DFN Jacobian mode: %s\n', mode);
    tic;

    sim = sim_profile_LGM50('DFN', init_soc, profile, ...
        'DeltaT', delta_t, ...
        'Tamb', t_amb, ...
        'Nr', Nr, 'Nxn', Nxn, 'Nxs', Nxs, 'Nxp', Nxp, ...
        'DegradationModel', 'safari', ...
        'EnableLithiumPlating', true, ...
        'EnableAnodePorosityLoss', true, ...
        'DFNLplCouplingMode', 'spme_compatible', ...
        'DFNJacobianMode', mode, ...
        'SaveOutput', false, ...
        'SuppressSummary', true);

    elapsed = toc;

    out.(mode).sim = sim;
    out.(mode).elapsed_s = elapsed;

    fprintf('Finished %s in %.2f s. Final V = %.6f, SOC = %.6f\n', ...
        mode, elapsed, sim.voltage(end), sim.soc(end));
end

plot_jacobian_mode_comparison(out, modes);

end

function plot_jacobian_mode_comparison(out, modes)

sim_fd = out.finite_difference.sim;
sim_pa = out.padded_analytic.sim;

t_fd = sim_fd.time(:);
t_pa = sim_pa.time(:);

n = min(numel(t_fd), numel(t_pa));
t = t_fd(1:n) / 3600;

v_fd = sim_fd.voltage(1:n);
v_pa = sim_pa.voltage(1:n);

soc_fd = sim_fd.soc(1:n);
soc_pa = sim_pa.soc(1:n);

d_fd = sim_fd.dfn_outputs;
d_pa = sim_pa.dfn_outputs;

j_fd = d_fd.j_lpl_sep_n(1:n);
j_pa = d_pa.j_lpl_sep_n(1:n);

li_fd = d_fd.delta_li_sep_n(1:n);
li_pa = d_pa.delta_li_sep_n(1:n);

eps_fd = d_fd.epsilon_e_sep_n(1:n);
eps_pa = d_pa.epsilon_e_sep_n(1:n);

ldp_fd = d_fd.LDP_sep_n(1:n);
ldp_pa = d_pa.LDP_sep_n(1:n);

fprintf('\nRuntime comparison:\n');
fprintf('finite_difference : %.3f s\n', out.finite_difference.elapsed_s);
fprintf('padded_analytic   : %.3f s\n', out.padded_analytic.elapsed_s);
fprintf('speedup           : %.2fx\n', out.finite_difference.elapsed_s / out.padded_analytic.elapsed_s);

fprintf('\nMax absolute trajectory differences:\n');
fprintf('Voltage      : %.3e V\n', max(abs(v_fd(:) - v_pa(:))));
fprintf('SOC          : %.3e\n', max(abs(soc_fd(:) - soc_pa(:))));
fprintf('LDP_sep      : %.3e V\n', max(abs(ldp_fd(:) - ldp_pa(:))));
fprintf('j_lpl_sep    : %.3e mol/m2/s\n', max(abs(j_fd(:) - j_pa(:))));
fprintf('delta_li_sep : %.3e m\n', max(abs(li_fd(:) - li_pa(:))));
fprintf('epsilon_sep  : %.3e\n', max(abs(eps_fd(:) - eps_pa(:))));

figure('Name','DFN Jacobian mode speed/accuracy comparison','Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on; box on;
bar(categorical(modes), [out.finite_difference.elapsed_s, out.padded_analytic.elapsed_s]);
ylabel('Runtime [s]');
title('Runtime');

nexttile; hold on; grid on; box on;
plot(t, v_fd, 'DisplayName','finite_difference');
plot(t, v_pa, '--', 'DisplayName','padded_analytic');
ylabel('Voltage [V]');
legend('Interpreter','none');

nexttile; hold on; grid on; box on;
plot(t, ldp_fd, 'DisplayName','finite_difference');
plot(t, ldp_pa, '--', 'DisplayName','padded_analytic');
ylabel('LDP_{sep} [V]');

nexttile; hold on; grid on; box on;
semilogy(t, max(abs(j_fd), realmin), 'DisplayName','finite_difference');
semilogy(t, max(abs(j_pa), realmin), '--', 'DisplayName','padded_analytic');
ylabel('|j_{lpl,sep}|');

nexttile; hold on; grid on; box on;
plot(t, li_fd*1e9, 'DisplayName','finite_difference');
plot(t, li_pa*1e9, '--', 'DisplayName','padded_analytic');
ylabel('\delta_{Li,sep} [nm]');
xlabel('Time [h]');

nexttile; hold on; grid on; box on;
plot(t, eps_fd, 'DisplayName','finite_difference');
plot(t, eps_pa, '--', 'DisplayName','padded_analytic');
ylabel('\epsilon_{e,sep}');
xlabel('Time [h]');

end
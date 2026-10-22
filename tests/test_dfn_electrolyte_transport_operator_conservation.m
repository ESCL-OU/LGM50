function out = test_dfn_electrolyte_transport_operator_conservation()
% Algebraic conservation test for DFN electrolyte transport operator.
%
% Purpose:
%   Test whether the electrolyte transport matrices conserve total
%   electrolyte Li inventory before running time integration.
%
% If W*M_transport is nonzero, the transport/boundary operator can create
% or remove electrolyte Li even with no SEI/plating.

init_soc = 0.1;
t_amb = 25;
delta_t = 1;
disable_thermal = true;

Nr = 10;
Nxn = 10;
Nxs = 5;
Nxp = 10;

p = load_DFN_params(init_soc, t_amb, delta_t, disable_thermal, ...
    Nr, Nxn, Nxs, Nxp, 'none', false, false, ...
    'spme_compatible', 'padded_analytic');

Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nce = p.Nx - 3;

if isfield(p, 'Area_ref')
    A = p.Area_ref;
else
    A = p.Area;
end

dx_n = p.L_n / Nn;
dx_s = p.L_s / Ns;
dx_p = p.L_p / Np;

W = [ ...
    p.epsilon_e_n * A * dx_n * ones(Nn,1); ...
    p.epsilon_e_s * A * dx_s * ones(Ns,1); ...
    p.epsilon_e_p * A * dx_p * ones(Np,1)];

Wrow = W(:).';

% Build linearized transport operator at uniform c_e.
% In dae_dfn:
% c_dot = dD.*grad^2 + D.*(M3*c + M4*c_bcs) + source
% For uniform c_e, dD term should be zero.
% We test the conservative part: M3 + M4*C.

E_n = [speye(Nn), sparse(Nn,Ns+Np)];
E_s = [sparse(Ns,Nn), speye(Ns), sparse(Ns,Np)];
E_p = [sparse(Np,Nn+Ns), speye(Np)];

M_transport_n = p.ce.M3n * E_n + p.ce.M4n * p.ce.C(1:2,:);
M_transport_s = p.ce.M3s * E_s + p.ce.M4s * p.ce.C(2:3,:);
M_transport_p = p.ce.M3p * E_p + p.ce.M4p * p.ce.C(3:4,:);

M_transport = [M_transport_n; M_transport_s; M_transport_p];

WM = Wrow * M_transport;

fprintf('\n=== DFN electrolyte transport conservation test ===\n');
fprintf('Nce = %d\n', Nce);
fprintf('norm(W*M_transport, inf) = %.6e\n', norm(WM, inf));
fprintf('norm(W*M_transport, 2)   = %.6e\n', norm(WM, 2));
fprintf('max abs W*M_transport    = %.6e\n', max(abs(WM)));
fprintf('sum W*M_transport        = %.6e\n', sum(WM));

% Test diffusion operator response to constant concentration.
c0 = p.c_e * ones(Nce,1);
transport_at_uniform = M_transport * c0;
int_transport_uniform = Wrow * transport_at_uniform;

fprintf('\nUniform c_e test:\n');
fprintf('max abs transport_at_uniform = %.6e mol/m^3/s\n', max(abs(transport_at_uniform)));
fprintf('integral transport uniform   = %.6e mol/s\n', int_transport_uniform);

% Test source conservation signs using actual initial z.
z0 = p.z0;
[~,~,~,~,~,jn0,jp0] = parse_z(z0,1,p);

source_n = p.ce.M5n * jn0;
source_p = p.ce.M5p * jp0;
source_full = [source_n; zeros(Ns,1); source_p];

int_source_n = Wrow(1:Nn) * source_n;
int_source_p = Wrow(Nn+Ns+1:end) * source_p;
int_source_total = Wrow * source_full;

fprintf('\nInitial source test:\n');
fprintf('integral source anode   = %.6e mol/s\n', int_source_n);
fprintf('integral source cathode = %.6e mol/s\n', int_source_p);
fprintf('integral source total   = %.6e mol/s\n', int_source_total);

% Finite perturbation tests
rng(1);
n_tests = 5;
fprintf('\nRandom c_e transport integral tests:\n');
random_integrals = zeros(n_tests,1);

for k = 1:n_tests
    c_rand = p.c_e * (1 + 0.2*randn(Nce,1));
    c_rand = max(c_rand, 1);
    val = Wrow * (M_transport * c_rand);
    random_integrals(k) = val;
    fprintf('test %d: W*M*c = %.6e mol/s\n', k, val);
end

% Localize offending columns
[abs_sorted, idx_sorted] = sort(abs(WM(:)), 'descend');

fprintf('\nLargest offending columns in W*M_transport:\n');
for k = 1:min(10,numel(idx_sorted))
    idx = idx_sorted(k);
    fprintf('col %3d | W*M = %+ .6e\n', idx, WM(idx));
end

% Package output
out = struct();
out.p = p;
out.W = W;
out.M_transport = M_transport;
out.WM = WM;
out.transport_at_uniform = transport_at_uniform;
out.int_transport_uniform = int_transport_uniform;
out.random_integrals = random_integrals;
out.source_initial = struct( ...
    'source_n', int_source_n, ...
    'source_p', int_source_p, ...
    'source_total', int_source_total);

plot_operator_conservation(out);

end

function plot_operator_conservation(out)

WM = out.WM(:);
transport_uniform = out.transport_at_uniform(:);

figure('Name','DFN electrolyte operator conservation test','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on; box on;
stem(WM, 'filled');
xlabel('c_e column index');
ylabel('W M_{transport}');
title('Weighted column sum of transport operator');

nexttile; hold on; grid on; box on;
stem(transport_uniform, 'filled');
xlabel('c_e node');
ylabel('transport(c_e = constant)');
title('Transport operator at uniform c_e');

nexttile; hold on; grid on; box on;
bar(out.random_integrals);
xlabel('random test');
ylabel('W M c [mol/s]');
title('Integrated transport for random c_e');

nexttile; hold on; grid on; box on;
bar([out.source_initial.source_n, out.source_initial.source_p, out.source_initial.source_total]);
set(gca,'XTickLabel',{'anode','cathode','total'});
ylabel('mol/s');
title('Initial integrated source');

end
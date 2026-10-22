function out = plot_jlpl_vs_eta_lpl_from_dfn(sim, varargin)
%PLOT_JLPL_VS_ETA_LPL_FROM_DFN
% Visualize Li plating flux/current density versus plating potential.
%
% Usage:
%   out = plot_jlpl_vs_eta_lpl_from_dfn(sim)
%   out = plot_jlpl_vs_eta_lpl_from_dfn(out.spme_compatible.sim)
%
% This uses the parameters in sim.p and overlays the actual DFN
% separator-side trajectory if sim.dfn_outputs exists.
%
% Sign convention:
%   j_lpl < 0 means cathodic Li plating flux.
%   eta_lpl is LDP / plating potential [V].
%
% Kinetic form:
%   j_lpl = -i0_lpl/F * exp((-alpha_lpl*F/(R*T))*eta_lpl)

parser = inputParser;
addParameter(parser, 'EtaMin', -0.20);
addParameter(parser, 'EtaMax',  0.20);
addParameter(parser, 'NEta', 2001);
addParameter(parser, 'UseAbsFlux', true);
addParameter(parser, 'ShowThresholdGate', true);
addParameter(parser, 'TitleSuffix', '');
parse(parser, varargin{:});
opt = parser.Results;

p = sim.p;

% Required params
i0_lpl = p.i0_lpl;      % [A/m^2]
alpha  = p.alph_lpl;
F      = p.Faraday;
R      = p.R;

if isfield(p, 'T_amb') && ~isempty(p.T_amb)
    T = p.T_amb;
else
    T = 298.15;
end

eta = linspace(opt.EtaMin, opt.EtaMax, opt.NEta).';  % [V]

% Irreversible Tafel plating flux [mol/m^2/s]
j_tafel = -(i0_lpl/F) .* exp((-alpha*F/(R*T)) .* eta);

% Threshold-gated version: zero unless eta < 0
j_gate = j_tafel;
j_gate(eta >= 0) = 0;

% Convert to current density [A/m^2]
i_tafel = F*j_tafel;
i_gate  = F*j_gate;

if opt.UseAbsFlux
    y_tafel = abs(j_tafel);
    y_gate  = abs(j_gate);
    y_label = '|j_{lpl}| [mol m^{-2} s^{-1}]';
else
    y_tafel = j_tafel;
    y_gate  = j_gate;
    y_label = 'j_{lpl} [mol m^{-2} s^{-1}]';
end

% Actual simulated trajectory, if available
has_dfn = isfield(sim, 'dfn_outputs') && ...
          isfield(sim.dfn_outputs, 'LDP_sep_n') && ...
          isfield(sim.dfn_outputs, 'j_lpl_sep_n');

eta_sim = [];
j_sim = [];
if has_dfn
    eta_sim = sim.dfn_outputs.LDP_sep_n(:);
    j_sim   = sim.dfn_outputs.j_lpl_sep_n(:);
end

if opt.UseAbsFlux && ~isempty(j_sim)
    j_sim_plot = abs(j_sim);
else
    j_sim_plot = j_sim;
end

% Print useful values around eta = 0
eta_probe = [-0.100; -0.050; -0.020; -0.010; -0.005; 0; 0.005; 0.010; 0.020; 0.050; 0.100];
j_probe = -(i0_lpl/F) .* exp((-alpha*F/(R*T)) .* eta_probe);

fprintf('\nLi plating Tafel curve using current DFN parameters:\n');
fprintf('  i0_lpl = %.6e A/m^2\n', i0_lpl);
fprintf('  alpha  = %.6g\n', alpha);
fprintf('  T      = %.6g K\n', T);
fprintf('  i0/F   = %.6e mol/m^2/s\n\n', i0_lpl/F);
fprintf(' eta_lpl [V]       j_lpl [mol/m^2/s]       |j_lpl|\n');
fprintf('----------------------------------------------------------\n');
for k = 1:numel(eta_probe)
    fprintf('% .6f        % .6e          %.6e\n', eta_probe(k), j_probe(k), abs(j_probe(k)));
end

if has_dfn
    fprintf('\nActual DFN separator-side trajectory:\n');
    fprintf('  eta_lpl range: [%.6e, %.6e] V\n', min(eta_sim), max(eta_sim));
    fprintf('  j_lpl range:   [%.6e, %.6e] mol/m^2/s\n', min(j_sim), max(j_sim));
    fprintf('  max |j_lpl|:   %.6e mol/m^2/s\n', max(abs(j_sim)));
end

% -------------------------
% Figure 1: wide eta range
% -------------------------
figure('Name', 'j_lpl vs eta_lpl', 'Color', 'w');
hold on; box on; grid on;

semilogy(eta, y_tafel, 'LineWidth', 1.8, 'DisplayName', 'Tafel, always active');

if opt.ShowThresholdGate
    semilogy(eta, max(y_gate, realmin), '--', 'LineWidth', 1.8, ...
        'DisplayName', 'Threshold-gated, eta < 0 only');
end

if has_dfn
    semilogy(eta_sim, max(j_sim_plot, realmin), '.', 'MarkerSize', 8, ...
        'DisplayName', 'DFN separator-side trajectory');
end

xline(0, 'k:', 'LineWidth', 1.2, 'DisplayName', '\eta_{lpl}=0');
xlabel('\eta_{lpl} / LDP [V]');
ylabel(y_label);

title_str = sprintf('Li plating flux vs plating potential%s', opt.TitleSuffix);
title(title_str, 'Interpreter', 'none');
legend('Location', 'best');

% -------------------------
% Figure 2: zoom around zero
% -------------------------
figure('Name', 'j_lpl vs eta_lpl near zero', 'Color', 'w');
hold on; box on; grid on;

idx_zoom = eta >= -0.03 & eta <= 0.03;

semilogy(eta(idx_zoom), y_tafel(idx_zoom), 'LineWidth', 1.8, ...
    'DisplayName', 'Tafel, always active');

if opt.ShowThresholdGate
    semilogy(eta(idx_zoom), max(y_gate(idx_zoom), realmin), '--', 'LineWidth', 1.8, ...
        'DisplayName', 'Threshold-gated, eta < 0 only');
end

if has_dfn
    idx_sim_zoom = eta_sim >= -0.03 & eta_sim <= 0.03;
    semilogy(eta_sim(idx_sim_zoom), max(j_sim_plot(idx_sim_zoom), realmin), '.', ...
        'MarkerSize', 8, 'DisplayName', 'DFN separator-side trajectory');
end

xline(0, 'k:', 'LineWidth', 1.2, 'DisplayName', '\eta_{lpl}=0');
xlabel('\eta_{lpl} / LDP [V]');
ylabel(y_label);
title('Zoom near \eta_{lpl}=0');
legend('Location', 'best');

% -------------------------
% Figure 3: eta and j over time
% -------------------------
if has_dfn && isfield(sim, 'time')
    t_hr = sim.time(:)/3600;
    n = min([numel(t_hr), numel(eta_sim), numel(j_sim)]);

    figure('Name', 'DFN separator-side eta_lpl and j_lpl over time', 'Color', 'w');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    nexttile;
    plot(t_hr(1:n), eta_sim(1:n), 'LineWidth', 1.4);
    hold on; grid on; box on;
    yline(0, 'k:', 'LineWidth', 1.2);
    ylabel('\eta_{lpl,sep} [V]');
    title('Separator-side plating potential over time');

    nexttile;
    if opt.UseAbsFlux
        semilogy(t_hr(1:n), max(abs(j_sim(1:n)), realmin), 'LineWidth', 1.4);
        ylabel('|j_{lpl,sep}| [mol m^{-2} s^{-1}]');
    else
        plot(t_hr(1:n), j_sim(1:n), 'LineWidth', 1.4);
        ylabel('j_{lpl,sep} [mol m^{-2} s^{-1}]');
    end
    grid on; box on;
    xlabel('Time [h]');
    title('Separator-side plating flux over time');
end

out = struct();
out.eta = eta;
out.j_tafel = j_tafel;
out.j_gate = j_gate;
out.i_tafel = i_tafel;
out.i_gate = i_gate;
out.eta_probe = eta_probe;
out.j_probe = j_probe;
out.has_dfn = has_dfn;
out.eta_sim = eta_sim;
out.j_sim = j_sim;

end
function out = plot_dfn_degradation_spatial_from_data(data, varargin)
%PLOT_DFN_DEGRADATION_SPATIAL_FROM_DATA
% Plot compact DFN spatial checkpoint data saved by
% test_dfn_degradation_spatial_checkpoints.m.
%
% Usage:
%   plot_dfn_degradation_spatial_from_data('results/DFN/dfn_degradation_spatial_checkpoints_yyyymmdd_HHMMSS.mat')
%   plot_dfn_degradation_spatial_from_data(checkpoint)
%
% X-axis is the DFN through-cell coordinate from negative current collector
% to positive current collector. Anode-only variables are plotted as NaN
% outside the anode.

    parser = inputParser;
    addParameter(parser, 'MaxCurves', 100);
    addParameter(parser, 'SaveFigures', true);
    addParameter(parser, 'FigureDir', '');
    addParameter(parser, 'UseAbsoluteFlux', true);
    parse(parser, varargin{:});
    opt = parser.Results;

    [checkpoint, source_label] = load_checkpoint_local(data);
    snapshots = checkpoint.snapshots;
    assert(~isempty(snapshots), 'No snapshots found.');

    idx = select_snapshot_indices_local(numel(snapshots), opt.MaxCurves);
    labels = cell(numel(idx), 1);
    for ii = 1:numel(idx)
        s = snapshots(idx(ii));
        labels{ii} = sprintf('cyc %d, seg %d', s.capacity_cycle, s.segment_idx);
    end

    if isempty(opt.FigureDir)
        if isfield(checkpoint, 'settings') && isfield(checkpoint.settings, 'run_timestamp')
            fig_dir = fullfile(pwd, 'results', 'DFN', ['figures_' checkpoint.settings.run_timestamp]);
        else
            fig_dir = fullfile(pwd, 'results', 'DFN', ['figures_' datestr(now, 'yyyymmdd_HHMMSS')]);
        end
    else
        fig_dir = opt.FigureDir;
    end
    if opt.SaveFigures && ~isfolder(fig_dir)
        mkdir(fig_dir);
    end

    % Figure 1: spatial electrochemical/degradation profiles
    fig1 = figure('Name', ['DFN spatial degradation checkpoints: ' source_label], 'Color', 'w');
    tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    plot_panel_local(1, snapshots, idx, labels, 'c_e', 'c_e [mol/m^3]', false, opt.UseAbsoluteFlux);
    plot_panel_local(2, snapshots, idx, labels, 'LDP', 'Plating potential / LDP [V]', false, opt.UseAbsoluteFlux);
    plot_panel_local(3, snapshots, idx, labels, 'delta_sei', '	heta_{SEI} or \delta_{SEI} [nm]', true, opt.UseAbsoluteFlux);
    plot_panel_local(4, snapshots, idx, labels, 'epsilon_e', '\epsilon_e [-]', false, opt.UseAbsoluteFlux);
    plot_panel_local(5, snapshots, idx, labels, 'delta_li', '	heta_{Li} or \delta_{Li} [nm]', true, opt.UseAbsoluteFlux);
    plot_panel_local(6, snapshots, idx, labels, 'j_lpl', 'j_{lpl} [mol m^{-2} s^{-1}]', false, opt.UseAbsoluteFlux);

    title(tl, 'DFN spatial profiles at saved checkpoints', 'Interpreter', 'none');

    % Figure 2: separator-side summaries over checkpoints
    fig2 = figure('Name', ['DFN separator-side checkpoint summary: ' source_label], 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    cycle = arrayfun(@(s) s.capacity_cycle, snapshots).';
    sep_LDP = arrayfun(@(s) s.separator_side.LDP, snapshots).';
    sep_j = arrayfun(@(s) s.separator_side.j_lpl, snapshots).';
    sep_dli = arrayfun(@(s) s.separator_side.delta_li, snapshots).';
    sep_eps = arrayfun(@(s) s.separator_side.epsilon_e, snapshots).';

    nexttile; plot(cycle, sep_LDP, 'o-', 'LineWidth', 1.1); grid on; xlabel('Capacity cycle'); ylabel('LDP_{sep} [V]');
    nexttile; plot(cycle, sep_j, 'o-', 'LineWidth', 1.1); grid on; xlabel('Capacity cycle'); ylabel('j_{lpl,sep} [mol m^{-2} s^{-1}]');
    nexttile; plot(cycle, sep_dli*1e9, 'o-', 'LineWidth', 1.1); grid on; xlabel('Capacity cycle'); ylabel('\delta_{Li,sep} [nm]');
    nexttile; plot(cycle, sep_eps, 'o-', 'LineWidth', 1.1); grid on; xlabel('Capacity cycle'); ylabel('\epsilon_{e,sep} [-]');

    if opt.SaveFigures
        savefig(fig1, fullfile(fig_dir, 'dfn_spatial_profiles.fig'));
        saveas(fig1, fullfile(fig_dir, 'dfn_spatial_profiles.png'));
        savefig(fig2, fullfile(fig_dir, 'dfn_separator_summary.fig'));
        saveas(fig2, fullfile(fig_dir, 'dfn_separator_summary.png'));
    end

    out = struct();
    out.checkpoint = checkpoint;
    out.source_label = source_label;
    out.selected_indices = idx;
    out.figure_dir = fig_dir;
end

function plot_panel_local(tile_idx, snapshots, idx, labels, field_name, y_label, scale_nm, use_abs_flux)
    nexttile(tile_idx); hold on; box on; grid on;
    for ii = 1:numel(idx)
        s = snapshots(idx(ii));
        y = s.(field_name);
        if scale_nm
            y = y * 1e9;
        end
        if strcmp(field_name, 'j_lpl') && use_abs_flux
            y = abs(y);
        end
        plot(s.x_norm, y, 'LineWidth', 1.1, 'DisplayName', labels{ii});
    end
    add_region_boundaries_local(snapshots(idx(end)));
    xlabel('Normalized through-cell position, x/L_{cell}');
    ylabel(y_label);
    title(field_name, 'Interpreter', 'none');
    if tile_idx == 1
        legend('Location', 'bestoutside', 'Interpreter', 'none');
    end
end

function add_region_boundaries_local(s)
    b = s.x_boundaries_norm;
    yl = ylim;
    xline(b(2), 'k--', 'HandleVisibility', 'off');
    xline(b(3), 'k--', 'HandleVisibility', 'off');
    ylim(yl);
end

function idx = select_snapshot_indices_local(n, max_curves)
    if n <= max_curves
        idx = 1:n;
    else
        idx = unique(round(linspace(1, n, max_curves)));
    end
end

function [checkpoint, source_label] = load_checkpoint_local(data)
    if ischar(data) || isstring(data)
        source_label = char(data);
        loaded = load(source_label);
        if isfield(loaded, 'checkpoint')
            checkpoint = loaded.checkpoint;
        elseif isfield(loaded, 'out')
            checkpoint = loaded.out;
        else
            error('No checkpoint variable found in %s.', source_label);
        end
    elseif isstruct(data) && isfield(data, 'snapshots')
        checkpoint = data;
        source_label = 'workspace checkpoint';
    elseif isstruct(data) && isfield(data, 'checkpoint')
        checkpoint = data.checkpoint;
        source_label = 'workspace checkpoint wrapper';
    else
        error('Input must be a checkpoint struct or path to .mat file.');
    end
end

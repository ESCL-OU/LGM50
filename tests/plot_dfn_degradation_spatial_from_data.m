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

    has_tags = isfield(snapshots, 'tag') && all(arrayfun(@(s) isfield(s, 'tag') && ~isempty(s.tag), snapshots));

    if has_tags
        idx_by_tag = group_snapshot_indices_by_tag_local(snapshots, opt.MaxCurves);
    else
        idx = select_snapshot_indices_local(numel(snapshots), opt.MaxCurves);
        idx_by_tag = struct();
        idx_by_tag.all = idx;
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
    fields_to_plot = { ...
        'c_e',       'c_e [mol/m^3]', false; ...
        'LDP',       'Plating potential / LDP [V]', false; ...
        'delta_sei', '\delta_{SEI} [nm]', true; ...
        'epsilon_e', '\epsilon_e [-]', false; ...
        'delta_li',  '\delta_{Li} [nm]', true; ...
        'j_lpl',     'j_{lpl} [mol m^{-2} s^{-1}]', false};
    
    tag_order = { ...
        'p0p2C_charge',    '0p2C_charge'; ...
        'p0p2C_discharge', '0p2C_discharge'; ...
        'p1C_charge',      '1C_charge'; ...
        'p1C_discharge',   '1C_discharge'};

    for f = 1:size(fields_to_plot,1)
        field_name = fields_to_plot{f,1};
        y_label = fields_to_plot{f,2};
        scale_nm = fields_to_plot{f,3};
    
        fig = figure('Name', ['DFN spatial ' field_name ': ' source_label], 'Color', 'w');
        tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
        for t = 1:size(tag_order,1)
            tag_key = tag_order{t,1};
            tag = tag_order{t,2};
            nexttile; hold on; box on; grid on;
    
            if isfield(idx_by_tag, tag_key)
                idx_this = idx_by_tag.(tag_key);
            else
                idx_this = [];
            end
    
            if isempty(idx_this)
                title(strrep(tag, '_', '\_'));
                xlabel('Normalized through-cell position, x/L_{cell}');
                ylabel(y_label);
                text(0.5, 0.5, 'No checkpoint', 'HorizontalAlignment', 'center');
                continue;
            end
    
            labels = make_snapshot_labels_local(snapshots, idx_this);
    
            for ii = 1:numel(idx_this)
                s = snapshots(idx_this(ii));
                y = s.(field_name);
    
                if scale_nm
                    y = y * 1e9;
                end
                if strcmp(field_name, 'j_lpl') && opt.UseAbsoluteFlux
                    y = abs(y);
                end
    
                plot(s.x_norm, y, 'LineWidth', 1.1, 'DisplayName', labels{ii});
            end
    
            add_region_boundaries_local(snapshots(idx_this(end)));
            xlabel('Normalized through-cell position, x/L_{cell}');
            ylabel(y_label);
            title(strrep(tag, '_', '\_'), 'Interpreter', 'tex');
    
            if numel(idx_this) <= 12
                legend('Location', 'best', 'Interpreter', 'none');
            end
        end
    
        title(tl, ['DFN spatial profiles: ' field_name], 'Interpreter', 'none');
    
        if opt.SaveFigures
            safe_name = regexprep(field_name, '[^A-Za-z0-9_]', '_');
            savefig(fig, fullfile(fig_dir, ['dfn_spatial_' safe_name '_by_tag.fig']));
            saveas(fig, fullfile(fig_dir, ['dfn_spatial_' safe_name '_by_tag.png']));
        end
    end
    % Figure 2: separator-side summaries over checkpoints
    fig2 = figure('Name', ['DFN separator-side checkpoint summary: ' source_label], 'Color', 'w');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    summary_fields = { ...
        'LDP',       'LDP_{sep} [V]'; ...
        'j_lpl',    'j_{lpl,sep} [mol m^{-2} s^{-1}]'; ...
        'delta_li', '\delta_{Li,sep} [nm]'; ...
        'epsilon_e','\epsilon_{e,sep} [-]'};
    
    for sf = 1:size(summary_fields,1)
        metric = summary_fields{sf,1};
        ylab = summary_fields{sf,2};
    
        nexttile; hold on; box on; grid on;
    
        for t = 1:numel(tag_order)
            tag = tag_order{t};
    
            if ~isfield(idx_by_tag, tag_key)
                continue;
            end
    
            idx_this = idx_by_tag.(tag_key);
            if isempty(idx_this)
                continue;
            end
    
            x_cycle = arrayfun(@(s) s.capacity_cycle, snapshots(idx_this)).';
    
            switch metric
                case 'LDP'
                    y = arrayfun(@(s) s.separator_side.LDP, snapshots(idx_this)).';
                case 'j_lpl'
                    y = arrayfun(@(s) s.separator_side.j_lpl, snapshots(idx_this)).';
                    if opt.UseAbsoluteFlux
                        y = abs(y);
                    end
                case 'delta_li'
                    y = arrayfun(@(s) s.separator_side.delta_li, snapshots(idx_this)).' * 1e9;
                case 'epsilon_e'
                    y = arrayfun(@(s) s.separator_side.epsilon_e, snapshots(idx_this)).';
            end
    
            plot(x_cycle, y, 'o-', 'LineWidth', 1.1, ...
                'DisplayName', strrep(tag, '_', ' '));
        end
    
        xlabel('Capacity cycle');
        ylabel(ylab);
        title(metric, 'Interpreter', 'none');
        legend('Location', 'best', 'Interpreter', 'none');
    end
    
    if opt.SaveFigures
        savefig(fig2, fullfile(fig_dir, 'dfn_separator_summary_by_tag.fig'));
        saveas(fig2, fullfile(fig_dir, 'dfn_separator_summary_by_tag.png'));
    end

    out = struct();
    out.checkpoint = checkpoint;
    out.source_label = source_label;
    if exist('idx', 'var')
        out.selected_indices = idx;
    else
        out.selected_indices_by_tag = idx_by_tag;
    end
    out.figure_dir = fig_dir;
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

function idx_by_tag = group_snapshot_indices_by_tag_local(snapshots, max_curves)
    tag_map = { ...
        'p0p2C_charge',    '0p2C_charge'; ...
        'p0p2C_discharge', '0p2C_discharge'; ...
        'p1C_charge',      '1C_charge'; ...
        'p1C_discharge',   '1C_discharge'};

    idx_by_tag = struct();

    for t = 1:size(tag_map,1)
        key = tag_map{t,1};
        tag = tag_map{t,2};

        idx = find(arrayfun(@(s) isfield(s, 'tag') && strcmp(s.tag, tag), snapshots));

        if numel(idx) > max_curves
            idx = idx(unique(round(linspace(1, numel(idx), max_curves))));
        end

        idx_by_tag.(key) = idx;
    end
end

function labels = make_snapshot_labels_local(snapshots, idx)
    labels = cell(numel(idx), 1);

    for ii = 1:numel(idx)
        s = snapshots(idx(ii));

        if isfield(s, 'sequence2_cycle_idx')
            labels{ii} = sprintf('cyc %d, seq2 %d, seg %d', ...
                s.capacity_cycle, s.sequence2_cycle_idx, s.segment_idx);
        elseif isfield(s, 'capacity_cycle')
            labels{ii} = sprintf('cyc %d, seg %d', ...
                s.capacity_cycle, s.segment_idx);
        else
            labels{ii} = sprintf('seg %d', s.segment_idx);
        end
    end
end
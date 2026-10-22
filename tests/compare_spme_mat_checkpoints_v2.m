function report = compare_spme_mat_checkpoints_v2(fileA, fileB, varargin)
%COMPARE_SPME_MAT_CHECKPOINTS Compare two SPMe degradation/checkpoint .mat files.
%
% Purpose
%   Use this when two runs look similar during degradation cycling (e.g., 1C),
%   but differ strongly at checkpoint / reference-capacity cycling (e.g., 0.2C
%   at cycle 200). The code recursively inspects both .mat files and reports:
%     1) top-level variables and fields that exist only in one file,
%     2) scalar/string/config differences,
%     3) numeric array differences with norm/max/relative metrics,
%     4) capacity-cycle-log differences if present,
%     5) terminal-state / checkpoint differences if present.
%
% Basic usage
%   report = compare_spme_mat_checkpoints('run_0625_cycle200.mat', ...
%                                        'run_0626_cycle200.mat');
%
% Suggested usage around your case
%   report = compare_spme_mat_checkpoints( ...
%       'spme_degradation_addons_checkpoint_safari_0625_repeat0002_segmentXXXX.mat', ...
%       'spme_degradation_addons_checkpoint_safari_0626_repeat0002_segmentXXXX.mat', ...
%       'CycleOfInterest', 200, ...
%       'TopN', 80, ...
%       'SavePrefix', 'diff_0625_vs_0626_cycle200');
%
% If you have many repeated checkpoint files, run compare_spme_checkpoint_dirs().
%
% Notes
%   - This is intentionally defensive: it does not assume an exact .mat schema.
%   - It can compare checkpoint structs, out structs, sim structs, or any saved variables.
%   - Large time-series arrays are summarized rather than printed in full.

p = inputParser;
p.addParameter('CycleOfInterest', 200, @(x) isempty(x) || isnumeric(x));
p.addParameter('TopN', 60, @(x) isnumeric(x) && isscalar(x));
p.addParameter('RelTol', 1e-9, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AbsTol', 1e-12, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MaxDepth', 12, @(x) isnumeric(x) && isscalar(x));
p.addParameter('SavePrefix', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

A = load(fileA);
B = load(fileB);

report = struct();
report.fileA = fileA;
report.fileB = fileB;
report.generated_at = datestr(now);
report.only_in_A = {};
report.only_in_B = {};
report.value_diffs = table();
report.numeric_diffs = table();
report.capacity_compare = table();
report.cycle_of_interest = opt.CycleOfInterest;

fprintf('\n=== Comparing SPMe MAT files ===\nA: %s\nB: %s\n', fileA, fileB);

% Recursive structural/numeric comparison
rowsNum = {};
rowsVal = {};
[rowsNum, rowsVal, onlyA, onlyB] = compare_any(A, B, '', rowsNum, rowsVal, {}, {}, 0, opt);
report.only_in_A = onlyA(:);
report.only_in_B = onlyB(:);
report.numeric_diffs = rows_to_numeric_table(rowsNum);
report.value_diffs = rows_to_value_table(rowsVal);

fprintf('\n--- Fields / variables only in A: %d, only in B: %d ---\n', numel(report.only_in_A), numel(report.only_in_B));
print_cell_list('Only in A', report.only_in_A, opt.TopN);
print_cell_list('Only in B', report.only_in_B, opt.TopN);

% Capacity log / cycle summary comparison if present
capA = find_capacity_table(A);
capB = find_capacity_table(B);
if ~isempty(capA) || ~isempty(capB)
    report.capacity_A = capA;
    report.capacity_B = capB;
    report.capacity_compare = compare_capacity_tables(capA, capB);
    fprintf('\n--- Capacity-cycle-log comparison ---\n');
    disp(limit_rows(report.capacity_compare, opt.TopN));

    if ~isempty(opt.CycleOfInterest) && ~isempty(report.capacity_compare)
        idx = report.capacity_compare.cycle == opt.CycleOfInterest;
        if any(idx)
            fprintf('\n*** Cycle %g capacity rows ***\n', opt.CycleOfInterest);
            disp(report.capacity_compare(idx, :));
        else
            fprintf('\nCycle %g was not found in the matched capacity table.\n', opt.CycleOfInterest);
        end
    end
end

% Terminal/checkpoint-focused summary
report.terminal_summary_A = summarize_checkpoint_like(A);
report.terminal_summary_B = summarize_checkpoint_like(B);
if ~isempty(report.terminal_summary_A) || ~isempty(report.terminal_summary_B)
    fprintf('\n--- Checkpoint / terminal summary A ---\n');
    disp(report.terminal_summary_A);
    fprintf('\n--- Checkpoint / terminal summary B ---\n');
    disp(report.terminal_summary_B);
end

% Print largest numeric differences
fprintf('\n--- Largest numeric differences ---\n');
if isempty(report.numeric_diffs)
    fprintf('No numeric differences found, or no comparable numeric arrays.\n');
else
    T = sortrows(report.numeric_diffs, {'rel_norm_diff','max_abs_diff'}, {'descend','descend'});
    disp(limit_rows(T, opt.TopN));
end

fprintf('\n--- Non-numeric/scalar/config differences ---\n');
if isempty(report.value_diffs)
    fprintf('No scalar/string/logical differences found.\n');
else
    disp(limit_rows(report.value_diffs, opt.TopN));
end

% Save report for later inspection
if strlength(string(opt.SavePrefix)) > 0
    matName = char(opt.SavePrefix + ".mat");
    csvNum = char(opt.SavePrefix + "_numeric_diffs.csv");
    csvVal = char(opt.SavePrefix + "_value_diffs.csv");
    save(matName, 'report');
    if ~isempty(report.numeric_diffs), writetable(report.numeric_diffs, csvNum); end
    if ~isempty(report.value_diffs), writetable(report.value_diffs, csvVal); end
    if ~isempty(report.capacity_compare)
        writetable(report.capacity_compare, char(opt.SavePrefix + "_capacity_compare.csv"));
    end
    fprintf('\nSaved report to: %s\n', matName);
end
end

function batch = compare_spme_checkpoint_dirs(dirA, dirB, varargin)
%COMPARE_SPME_CHECKPOINT_DIRS Pair checkpoint .mat files from two run folders.
%
% Usage
%   batch = compare_spme_checkpoint_dirs('results/SPMe/0625', 'results/SPMe/0626', ...
%       'Pattern', '*checkpoint*.mat', 'CycleOfInterest', 200);
%
% Pairing rule
%   First tries to pair by repeat####_segment##### in filename. If absent, pairs
%   sorted files by index.

p = inputParser;
p.addParameter('Pattern', '*checkpoint*.mat', @(x) ischar(x) || isstring(x));
p.addParameter('CycleOfInterest', 200, @(x) isempty(x) || isnumeric(x));
p.addParameter('TopN', 30, @(x) isnumeric(x) && isscalar(x));
p.addParameter('SavePrefix', 'batch_diff_0625_vs_0626', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

filesA = dir(fullfile(dirA, opt.Pattern));
filesB = dir(fullfile(dirB, opt.Pattern));
filesA = filesA(~[filesA.isdir]);
filesB = filesB(~[filesB.isdir]);

keysA = arrayfun(@(f) checkpoint_key(f.name), filesA, 'UniformOutput', false);
keysB = arrayfun(@(f) checkpoint_key(f.name), filesB, 'UniformOutput', false);

pairs = {};
usedB = false(numel(filesB), 1);
for i = 1:numel(filesA)
    j = find(strcmp(keysB, keysA{i}) & ~usedB, 1, 'first');
    if isempty(j) && i <= numel(filesB) && ~usedB(i)
        j = i;
    end
    if ~isempty(j)
        pairs(end+1, :) = {fullfile(filesA(i).folder, filesA(i).name), fullfile(filesB(j).folder, filesB(j).name), keysA{i}}; %#ok<AGROW>
        usedB(j) = true;
    end
end

batch = table('Size', [size(pairs,1), 8], ...
    'VariableTypes', {'string','string','string','double','double','double','double','double'}, ...
    'VariableNames', {'fileA','fileB','key','n_numeric_diffs','max_abs_diff','max_rel_norm_diff','cycle200_discharged_Ah_A','cycle200_discharged_Ah_B'});

for k = 1:size(pairs, 1)
    fprintf('\n\n######## Pair %d / %d: %s ########\n', k, size(pairs, 1), pairs{k,3});
    r = compare_spme_mat_checkpoints(pairs{k,1}, pairs{k,2}, ...
        'CycleOfInterest', opt.CycleOfInterest, 'TopN', opt.TopN);
    batch.fileA(k) = string(pairs{k,1});
    batch.fileB(k) = string(pairs{k,2});
    batch.key(k) = string(pairs{k,3});
    batch.n_numeric_diffs(k) = height(r.numeric_diffs);
    if ~isempty(r.numeric_diffs)
        batch.max_abs_diff(k) = max(r.numeric_diffs.max_abs_diff, [], 'omitnan');
        batch.max_rel_norm_diff(k) = max(r.numeric_diffs.rel_norm_diff, [], 'omitnan');
    end
    if ~isempty(r.capacity_compare) && ~isempty(opt.CycleOfInterest)
        idx = r.capacity_compare.cycle == opt.CycleOfInterest;
        if any(idx)
            if ismember('discharged_Ah_A', r.capacity_compare.Properties.VariableNames)
                batch.cycle200_discharged_Ah_A(k) = r.capacity_compare.discharged_Ah_A(find(idx,1));
            end
            if ismember('discharged_Ah_B', r.capacity_compare.Properties.VariableNames)
                batch.cycle200_discharged_Ah_B(k) = r.capacity_compare.discharged_Ah_B(find(idx,1));
            end
        end
    end
end

writetable(batch, char(opt.SavePrefix + ".csv"));
save(char(opt.SavePrefix + ".mat"), 'batch');
fprintf('\nSaved batch summary to %s.csv and %s.mat\n', opt.SavePrefix, opt.SavePrefix);
end

function key = checkpoint_key(name)
tok = regexp(name, '(repeat\d+_segment\d+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(name, '(cycle\d+)', 'tokens', 'once');
end
if isempty(tok)
    key = name;
else
    key = tok{1};
end
end

function [rowsNum, rowsVal, onlyA, onlyB] = compare_any(a, b, path, rowsNum, rowsVal, onlyA, onlyB, depth, opt)
if depth > opt.MaxDepth
    return;
end

if isstruct(a) && isstruct(b)
    if numel(a) ~= numel(b)
        rowsVal(end+1,:) = {path, class(a), class(b), size_string(a), size_string(b), 'struct array size differs'};
    end
    if numel(a) ~= 1 || numel(b) ~= 1
        n = min(numel(a), numel(b));
        for i = 1:min(n, 5)
            [rowsNum, rowsVal, onlyA, onlyB] = compare_any(a(i), b(i), sprintf('%s(%d)', path, i), rowsNum, rowsVal, onlyA, onlyB, depth+1, opt);
        end
        return;
    end
    fa = fieldnames(a); fb = fieldnames(b);
    onlyA = [onlyA; strcat(path, '.', setdiff(fa, fb))]; %#ok<AGROW>
    onlyB = [onlyB; strcat(path, '.', setdiff(fb, fa))]; %#ok<AGROW>
    common = intersect(fa, fb);
    for i = 1:numel(common)
        f = common{i};
        child = f;
        if ~isempty(path), child = [path '.' f]; end %#ok<AGROW>
        [rowsNum, rowsVal, onlyA, onlyB] = compare_any(a.(f), b.(f), child, rowsNum, rowsVal, onlyA, onlyB, depth+1, opt);
    end
    return;
end

if istable(a) && istable(b)
    [rowsNum, rowsVal, onlyA, onlyB] = compare_any(table2struct(a, 'ToScalar', true), table2struct(b, 'ToScalar', true), path, rowsNum, rowsVal, onlyA, onlyB, depth+1, opt);
    return;
end

if isnumeric(a) && isnumeric(b)
    rowsNum = compare_numeric(a, b, path, rowsNum, opt);
    return;
end

if islogical(a) && islogical(b)
    rowsNum = compare_numeric(double(a), double(b), path, rowsNum, opt);
    return;
end

if (ischar(a) || isstring(a)) && (ischar(b) || isstring(b))
    if ~isequal(string(a), string(b))
        rowsVal(end+1,:) = {path, class(a), class(b), char(string(a)), char(string(b)), 'string differs'}; %#ok<AGROW>
    end
    return;
end

if ~strcmp(class(a), class(b)) || ~isequal(size(a), size(b))
    rowsVal(end+1,:) = {path, class(a), class(b), size_string(a), size_string(b), 'class or size differs'}; %#ok<AGROW>
end
end

function rowsNum = compare_numeric(a, b, path, rowsNum, opt)
if ~isequal(size(a), size(b))
    rowsNum(end+1,:) = {path, class(a), size_string(a), size_string(b), NaN, NaN, NaN, NaN, NaN, 'numeric size differs'};
    return;
end
if isempty(a) && isempty(b), return; end
ad = double(a(:));
bd = double(b(:));

% Treat matching NaNs as equal, but count one-sided NaNs/Infs as differences.
sameNaN = isnan(ad) & isnan(bd);
idx = ~sameNaN;
if ~any(idx)
    return;
end

av = ad(idx);
bv = bd(idx);
d = av - bv;
absD = abs(d);

finiteD = absD(isfinite(absD));
if isempty(finiteD)
    maxAbs = double(any(absD ~= 0 | xor(isnan(av), isnan(bv)) | xor(isinf(av), isinf(bv))));
    meanAbs = maxAbs;
else
    maxAbs = max(finiteD);
    meanAbs = mean(finiteD);
end

finiteDiff = d(isfinite(d));
if isempty(finiteDiff)
    normDiff = maxAbs;
else
    normDiff = norm(finiteDiff);
end

finiteA = ad(isfinite(ad));
finiteB = bd(isfinite(bd));
if isempty(finiteA), normA = 0; else, normA = norm(finiteA); end
if isempty(finiteB), normB = 0; else, normB = norm(finiteB); end
relNorm = normDiff / max([normA, normB, eps]);

threshold = opt.AbsTol + opt.RelTol .* max(abs(av), abs(bv));
changedMask = absD > threshold;
changedMask = changedMask | xor(isnan(av), isnan(bv)) | xor(isinf(av), isinf(bv));
nDiff = sum(changedMask(:));

% Force scalar numeric metrics so MATLAB if/|| never sees vectors or [].
maxAbs = double(maxAbs(1));
meanAbs = double(meanAbs(1));
normDiff = double(normDiff(1));
relNorm = double(relNorm(1));
nDiff = double(nDiff(1));

absTolScalar = double(max(opt.AbsTol(:)));
relTolScalar = double(max(opt.RelTol(:)));
if any(nDiff(:) > 0) || any(maxAbs(:) > absTolScalar) || any(relNorm(:) > relTolScalar)
    rowsNum(end+1,:) = {path, class(a), size_string(a), size_string(b), double(nDiff(1)), double(maxAbs(1)), double(meanAbs(1)), double(normDiff(1)), double(relNorm(1)), ''}; %#ok<AGROW>
end
end

function T = rows_to_numeric_table(rows)
vars = {'path','class','size_A','size_B','n_changed','max_abs_diff','mean_abs_diff','norm_diff','rel_norm_diff','note'};
if isempty(rows)
    T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    T = cell2table(rows, 'VariableNames', vars);
end
end

function T = rows_to_value_table(rows)
vars = {'path','class_A','class_B','value_A','value_B','note'};
if isempty(rows)
    T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    T = cell2table(rows, 'VariableNames', vars);
end
end

function s = size_string(x)
sz = size(x);
s = strjoin(string(sz), 'x');
s = char(s);
end

function print_cell_list(titleStr, c, topN)
fprintf('%s:\n', titleStr);
if isempty(c)
    fprintf('  <none>\n');
else
    for i = 1:min(numel(c), topN)
        fprintf('  %s\n', c{i});
    end
    if numel(c) > topN, fprintf('  ... %d more\n', numel(c)-topN); end
end
end

function T = limit_rows(T, n)
if height(T) > n
    T = T(1:n, :);
end
end

function cap = find_capacity_table(S)
cap = table();
paths = {'capacity_cycle_log', 'sim.capacity_cycle_log', 'out.cycle_id'};
% Direct capacity_cycle_log struct array
for i = 1:numel(paths)
    val = get_by_path(S, paths{i});
    if ~isempty(val)
        cap = capacity_to_table(val, S);
        if ~isempty(cap), return; end
    end
end
% Recursive search for a field named capacity_cycle_log
val = recursive_find_field(S, 'capacity_cycle_log', 0, 8);
if ~isempty(val)
    cap = capacity_to_table(val, S);
end
end

function T = capacity_to_table(x, root)
T = table();
if isstruct(x) && all(isfield(x, {'cycle','charged_Ah','discharged_Ah'}))
    T = struct2table(x(:));
    keep = intersect({'cycle','charged_Ah','discharged_Ah','net_Ah'}, T.Properties.VariableNames, 'stable');
    T = T(:, keep);
elseif isnumeric(x)
    % Fallback: if only out.cycle_id was found, try sibling fields.
    cy = x(:);
    ch = get_by_path(root, 'out.charged_ah');
    dis = get_by_path(root, 'out.discharged_ah');
    net = get_by_path(root, 'out.net_ah');
    if ~isempty(ch) && ~isempty(dis)
        T = table(cy, ch(:), dis(:), 'VariableNames', {'cycle','charged_Ah','discharged_Ah'});
        if ~isempty(net), T.net_Ah = net(:); end
    end
end
end

function C = compare_capacity_tables(A, B)
C = table();
if isempty(A) || isempty(B), return; end
if ~ismember('cycle', A.Properties.VariableNames) || ~ismember('cycle', B.Properties.VariableNames), return; end
common = intersect(A.cycle, B.cycle);
C = table(common, 'VariableNames', {'cycle'});
vars = setdiff(intersect(A.Properties.VariableNames, B.Properties.VariableNames), {'cycle'}, 'stable');
for i = 1:numel(vars)
    v = vars{i};
    [~, ia] = ismember(common, A.cycle);
    [~, ib] = ismember(common, B.cycle);
    av = A.(v)(ia); bv = B.(v)(ib);
    C.([v '_A']) = av;
    C.([v '_B']) = bv;
    C.([v '_diff_B_minus_A']) = bv - av;
    C.([v '_rel_diff_B_minus_A']) = (bv - av) ./ max(abs(av), eps);
end
end

function summary = summarize_checkpoint_like(S)
summary = table();
rows = {};
paths = {'checkpoint.terminal_soc', 'checkpoint.terminal_voltage', 'checkpoint.next_segment_idx', ...
         'checkpoint.terminal_state.x', 'terminal_state.x', 'x', ...
         'out.cycle_id', 'out.discharged_ah', 'out.charged_ah', ...
         'sim.time', 'sim.current', 'sim.voltage', 'sim.soc'};
for i = 1:numel(paths)
    v = get_by_path(S, paths{i});
    if ~isempty(v)
        rows(end+1,:) = summarize_value(paths{i}, v); %#ok<AGROW>
    end
end
if ~isempty(rows)
    summary = cell2table(rows, 'VariableNames', {'path','class','size','first','last','min','max','mean'});
end
end

function row = summarize_value(path, v)
if isnumeric(v) || islogical(v)
    vv = double(v(:));
    row = {path, class(v), size_string(v), first_num(vv), last_num(vv), min(vv, [], 'omitnan'), max(vv, [], 'omitnan'), mean(vv, 'omitnan')};
elseif ischar(v) || isstring(v)
    row = {path, class(v), size_string(v), char(string(v)), char(string(v)), NaN, NaN, NaN};
else
    row = {path, class(v), size_string(v), NaN, NaN, NaN, NaN, NaN};
end
end

function y = first_num(v)
if isempty(v), y = NaN; else, y = v(1); end
end
function y = last_num(v)
if isempty(v), y = NaN; else, y = v(end); end
end

function v = get_by_path(S, path)
v = [];
parts = split(string(path), '.');
cur = S;
for i = 1:numel(parts)
    f = char(parts(i));
    if isstruct(cur) && isfield(cur, f)
        cur = cur.(f);
    else
        return;
    end
end
v = cur;
end

function val = recursive_find_field(S, target, depth, maxDepth)
val = [];
if depth > maxDepth || ~isstruct(S), return; end
if isfield(S, target)
    val = S.(target);
    return;
end
fn = fieldnames(S);
for i = 1:numel(fn)
    child = S.(fn{i});
    if isstruct(child)
        val = recursive_find_field(child, target, depth+1, maxDepth);
        if ~isempty(val), return; end
    end
end
end

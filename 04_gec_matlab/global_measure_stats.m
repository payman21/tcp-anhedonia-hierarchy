% File: analyze_pvalues_from_saved_results.m
% Purpose: compute between-group p-values from outputs saved by the original script.
% Inputs (existing files produced by the original code):
%   - outputs/results_Ceff_anhedonia_comparison.mat   [preferred]
%     or
%   - outputs/low_anhedonia/results_Ceff_low_anhedonia.mat
%   - outputs/high_anhedonia/results_Ceff_high_anhedonia.mat
% Outputs:
%   - prints p-values to console
%   - saves outputs/pvals_summary.mat and outputs/pvals_summary.csv

clear; clc;

%% Locate and load saved results (prefer combined file)
base_dir = 'outputs';
combo_fp = fullfile(base_dir, 'results_Ceff_anhedonia_comparison.mat');
low_fp   = fullfile(base_dir, 'low_anhedonia',  'results_Ceff_low_anhedonia.mat');
high_fp  = fullfile(base_dir, 'high_anhedonia', 'results_Ceff_high_anhedonia.mat');

if exist(combo_fp, 'file')
    S = load(combo_fp);
elseif exist(low_fp, 'file') && exist(high_fp, 'file')
    Slo = load(low_fp);
    Shi = load(high_fp);
    % harmonize into struct S (only required fields)
    S.trophiccoherence_LOW  = Slo.trophiccoherence_LOW;
    S.trophiccoherence_HIGH = Shi.trophiccoherence_HIGH;
    S.hierarchicallevels_LOW  = Slo.hierarchicallevels_LOW;
    S.hierarchicallevels_HIGH = Shi.hierarchicallevels_HIGH;
    S.fittFC_LOW   = Slo.fittFC_LOW;
    S.fittFC_HIGH  = Shi.fittFC_HIGH;
    S.fittCVtau_LOW  = Slo.fittCVtau_LOW;
    S.fittCVtau_HIGH = Shi.fittCVtau_HIGH;
else
    error('Could not find saved results. Expected %s or the two group files.', combo_fp);
end

%% Extract subject-level vectors
% Correlation-based fits
fc_low     = S.fittFC_LOW(:);
fc_high    = S.fittFC_HIGH(:);
cvtau_low  = S.fittCVtau_LOW(:);
cvtau_high = S.fittCVtau_HIGH(:);

% Trophic coherence (already compared in original code; recompute here)
tc_low  = S.trophiccoherence_LOW(:);
tc_high = S.trophiccoherence_HIGH(:);

% Mean hierarchical level per subject
mh_low  = mean(S.hierarchicallevels_LOW,  2); mh_low  = mh_low(:);
mh_high = mean(S.hierarchicallevels_HIGH, 2); mh_high = mh_high(:);

%% Nonparametric tests (Wilcoxon rank-sum / Mann–Whitney U)
[p_tc_np,  ~] = ranksum(tc_low,  tc_high);
[p_fc_np,  ~] = ranksum(fc_low,  fc_high);
[p_tau_np, ~] = ranksum(cvtau_low, cvtau_high);
[p_mh_np,  ~] = ranksum(mh_low,  mh_high);

p_np = [p_tc_np; p_fc_np; p_tau_np; p_mh_np];

%% Parametric tests (Welch's t)
% Fisher z-transform for correlations before parametric testing
z = @(r) atanh(max(min(r, 0.999999), -0.999999)); % clamp to avoid inf
[~, p_fc_w]  = ttest2(z(fc_low),   z(fc_high),   'Vartype','unequal');
[~, p_tau_w] = ttest2(z(cvtau_low),z(cvtau_high),'Vartype','unequal');
[~, p_mh_w]  = ttest2(mh_low, mh_high, 'Vartype','unequal');
% For completeness, trophic coherence with Welch as well
[~, p_tc_w]  = ttest2(tc_low, tc_high, 'Vartype','unequal');

p_w = [p_tc_w; p_fc_w; p_tau_w; p_mh_w];

%% Multiple-comparison control (Benjamini–Hochberg FDR; Bonferroni)
bh = @(p) benjamini_hochberg(p);
p_np_fdr = bh(p_np);
p_w_fdr  = bh(p_w);

m = numel(p_np);
p_np_bonf = min(1, p_np*m);
p_w_bonf  = min(1, p_w*m);

%% Assemble table
labels = {'TrophicCoherence','FC_fit_corr','COVtau_fit_corr','MeanHierLevel'}';
T = table(labels, ...
    p_np, p_np_fdr, p_np_bonf, ...
    p_w,  p_w_fdr,  p_w_bonf, ...
    'VariableNames', {'Measure','p_Wilcoxon','p_Wilcoxon_FDR','p_Wilcoxon_Bonf', ...
                                  'p_Welch','p_Welch_FDR','p_Welch_Bonf'});

%% Print concise report
disp('Between-group p-values (Low vs High):');
disp(T);

%% Save
if ~exist(base_dir,'dir'), mkdir(base_dir); end
save(fullfile(base_dir,'pvals_summary.mat'), 'T');
writetable(T, fullfile(base_dir,'pvals_summary.csv'));

%% ---- Helper: Benjamini–Hochberg FDR ----
function p_adj = benjamini_hochberg(p)
    [p_sorted, idx] = sort(p(:));
    m = numel(p_sorted);
    ranks = (1:m)';
    q = p_sorted .* m ./ ranks;
    % enforce monotonicity
    for i = m-1:-1:1
        q(i) = min(q(i), q(i+1));
    end
    p_adj = zeros(size(p_sorted));
    p_adj(:) = min(q, 1);
    % unsort
    p_adj(idx) = p_adj;
    p_adj = reshape(p_adj, size(p));
end

%%


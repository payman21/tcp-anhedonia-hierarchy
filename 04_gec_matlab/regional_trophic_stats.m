% % File: plot_pvalues.m
% % Purpose: make a significance plot from a vector of p-values only.
% % Inputs you must provide in the workspace OR replace below:
% %   p_vals : [N x 1] or [1 x N] vector of p-values for N tests
% % Optional:
% %   labels : {N x 1} cellstr of region names (used only for hover/print)
% %   outdir : output folder (default 'outputs')
% 
% function plot_pvalues(p_vals, labels, outdir)
%     if nargin < 2 || isempty(labels), labels = {}; end
%     if nargin < 3 || isempty(outdir), outdir = 'outputs'; end
%     if ~isvector(p_vals), error('p_vals must be a vector'); end
%     p_vals = p_vals(:);
%     N = numel(p_vals);
% 
%     % Clean/clip p-values to avoid Inf in -log10
%     eps_min = 1e-300;
%     p_clipped = max(min(p_vals, 1), eps_min);
% 
%     % Multiple-comparison thresholds
%     [p_fdr, thr_bh] = bh_fdr(p_clipped, 0.05);
%     sig_bh = p_fdr < 0.05;
% 
%     % Bonferroni
%     alpha = 0.05;
%     alpha_bonf = alpha / N;
% 
%     % Scatter plot of -log10(p)
%     x = (1:N)';
%     y = -log10(p_clipped);
% 
%     figure('Color','w','Position',[100 100 900 500]);
%     scatter(x, y, 26, 'k', 'filled'); hold on;
% 
%     % Threshold lines
%     y_nom  = -log10(alpha);
%     y_bonf = -log10(alpha_bonf);
%     y_bh   = -log10(thr_bh);
% 
%         % Add this after calculating the thresholds:
%     fprintf('Nominal p=0.05 threshold: y = %.2f\n', y_nom);
%     fprintf('Bonferroni threshold: y = %.2f\n', y_bonf); 
%     fprintf('BH FDR threshold: y = %.2f (raw p = %.2e)\n', y_bh, thr_bh);
% 
%     yline(y_nom,  '--r', 'LineWidth', 1);      % nominal 0.05
%     yline(y_bh,   '--b', 'LineWidth', 1);      % BH FDR 0.05 cutoff
%     yline(y_bonf, '--k', 'LineWidth', 1);      % Bonferroni
% 
%     % Highlight FDR-significant points
%     scatter(x(sig_bh), y(sig_bh), 34, 'b', 'filled', 'MarkerEdgeColor','w');
% 
%     % Axis labels, title, legend with 200% font size
%     ax = gca;
%     ax.FontSize = 20;  % tick labels
% 
%     xlabel('Brain Region Index','FontSize',20);
%     ylabel('-log_{10}(p)','FontSize',20);
%     title('P-value map with FDR and Bonferroni thresholds','FontSize',22);
%     xlim([1 232]);  % force x-axis to stop at 232
% 
%     lgd = legend({'tests', 'p=0.05', 'BH FDR q=0.05', 'Bonferroni'}, 'Location','best');
%     lgd.FontSize = 18;
% 
%     box off; grid on;
% 
%     % Optional: print list of significant indices/names
%     if any(sig_bh)
%         fprintf('FDR<0.05 indices: %s\n', mat2str(find(sig_bh)'));
%         if ~isempty(labels) && numel(labels)==N
%             disp('FDR<0.05 names:'); disp(labels(sig_bh));
%         end
%     else
%         fprintf('No tests survive FDR<0.05\n');
%     end
% 
%     % Save
%     if ~exist(outdir,'dir'), mkdir(outdir); end
%     savefig(fullfile(outdir,'pvalue_map.fig'));
%     print(fullfile(outdir,'pvalue_map.png'), '-dpng', '-r300');
% end
% 
% % ---------- Helpers ----------
% function [p_adj, thr] = bh_fdr(p, q)
%     % Benjamini–Hochberg adjustment (returns adjusted p-values and cutoff)
%     [ps, idx] = sort(p(:));
%     m = numel(ps);
%     ranks = (1:m)';
%     qline = (ranks/m) * q;
%     below = ps <= qline;
%     if any(below)
%         k = find(below, 1, 'last');
%         thr = ps(k);           % BH threshold (on raw p scale)
%     else
%         thr = 0;               % no discoveries
%     end
%     % adjusted p-values
%     adj = ps .* m ./ ranks;
%     adj = cummin(flipud(adj));
%     adj = flipud(adj);
%     adj = min(adj, 1);
%     p_adj = zeros(size(p));
%     p_adj(idx) = adj;
% end
% 
% 
% plot_pvalues(p_vals);  

% --- 0. LOAD DATA ---
Slo = load(fullfile('outputs', 'low_anhedonia', 'results_Ceff_low_anhedonia.mat'), ...
           'hierarchicallevels_LOW');
Shi = load(fullfile('outputs', 'high_anhedonia', 'results_Ceff_high_anhedonia.mat'), ...
           'hierarchicallevels_HIGH');
hierarchicallevels_LOW  = Slo.hierarchicallevels_LOW;
hierarchicallevels_HIGH = Shi.hierarchicallevels_HIGH;

% --- 1. CALCULATE P-VALUES PER REGION ---
% We will compare the hierarchical levels for each region between groups.
N_regions = size(hierarchicallevels_LOW, 2);
p_values_regions = zeros(1, N_regions);

for i = 1:N_regions
    % Compare the i-th column (region) of Low vs High
    p_values_regions(i) = ranksum(hierarchicallevels_LOW(:, i), ...
                                  hierarchicallevels_HIGH(:, i));
end

% Convert to -log10 scale for plotting
neg_log_p = -log10(p_values_regions);

% --- 2. CALCULATE THRESHOLDS ---
alpha = 0.05;

% A. Uncorrected Threshold (p = 0.05)
thresh_unc = -log10(alpha);

% B. Bonferroni Threshold (p = 0.05 / N)
thresh_bonf = -log10(alpha / N_regions);

% C. FDR Benjamini-Hochberg Threshold
% Sort p-values
[p_sorted, sort_idx] = sort(p_values_regions);
% Find largest k such that p(k) <= (k/N)*alpha
k_vec = 1:N_regions;
fdr_line = (k_vec / N_regions) * alpha;
is_sig = p_sorted <= fdr_line;

if any(is_sig)
    % The threshold is the p-value of the last significant test
    max_k = find(is_sig, 1, 'last');
    fdr_p_val = p_sorted(max_k);
    thresh_fdr = -log10(fdr_p_val);
else
    % If nothing is FDR significant, we can't draw a line, or we draw it at 0
    % For visualization, sometimes people plot the theoretical line, but a
    % horizontal cutoff is only valid if at least one value passes.
    thresh_fdr = NaN; 
end

% --- 3. GENERATE THE MANHATTAN PLOT ---
figure('Color', 'w', 'Position', [100, 100, 900, 500]); hold on;

% Plot the p-values (Black dots)
scatter(1:N_regions, neg_log_p, 25, 'k', 'filled', 'DisplayName', 'tests');

% Plot Uncorrected Threshold (Red dashed)
yline(thresh_unc, 'r--', 'p=0.05', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', 'DisplayName', 'p=0.05');

% Plot FDR Threshold (Blue dashed) - Only if valid
if ~isnan(thresh_fdr)
    yline(thresh_fdr, 'b--', 'BH FDR q=0.05', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', 'DisplayName', 'BH FDR q=0.05');
end

% Plot Bonferroni Threshold (Black dashed)
yline(thresh_bonf, 'k--', 'Bonferroni', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', 'DisplayName', 'Bonferroni');

% --- 4. FORMATTING ---
xlabel('Brain Region Index', 'FontSize', 14);
ylabel('-log_{10}(p)', 'FontSize', 14);
title('P-value map with FDR and Bonferroni thresholds', 'FontSize', 16, 'FontWeight', 'bold');
xlim([0, N_regions+5]);
ylim([0, max([max(neg_log_p), thresh_bonf]) * 1.1]); % Scale Y to fit data + bonferroni
grid on;
set(gca, 'FontSize', 12);
legend('Location', 'best');

hold off;
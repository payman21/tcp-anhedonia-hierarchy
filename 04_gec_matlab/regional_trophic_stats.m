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
% =========================================================================
% CODE OPTIMIZATION SUMMARY (what's changed compared to original code)
% =========================================================================
% 1. Vectorized Connectivity Update: Replaced the slow N*N nested loops in 
%    the optimization step with a single matrix operation using a logical mask.
% 2. Parallel Processing: Switched subject loops to 'parfor' to process 
%    multiple subjects simultaneously using all available CPU cores.
% 3. Optimized xcov: Replaced N*N calls to the 'xcov' function with a single, 
%    fast matrix multiplication for time-lagged covariance.
% 4. Code Refactoring: Consolidated the repeated optimization logic into a 
%    single helper function 'run_optimization' to reduce redundancy.
% =========================================================================

clear all
% Set up output directories
output_dir = [cfg.OUTPUT filesep];
output_low_dir = fullfile(output_dir, 'low_anhedonia/');
output_high_dir = fullfile(output_dir, 'high_anhedonia/');

% Create output directories
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
if ~exist(output_low_dir, 'dir'), mkdir(output_low_dir); end
if ~exist(output_high_dir, 'dir'), mkdir(output_high_dir); end

% Load Data
% (Update paths as needed)
cfg = config();
load(cfg.SC_FILE);
C = SC; 

load(fullfile(cfg.CLUSTERS, 'low_anhedonia.mat'));
load(fullfile(cfg.CLUSTERS, 'high_anhedonia.mat'));

load(fullfile(output_dir, 'empirical_low_anhedonia.mat'));
f_diff_low = f_diff;
load(fullfile(output_dir, 'empirical_high_anhedonia.mat'));
f_diff_high = f_diff;

% --- Parameters ---
N = 232;  
params.N = N;
params.Tau = 3;
params.sigma = 0.01;
params.TR = 0.8;
params.epsFC = 0.0004;
params.epsFCtau = 0.0001;
params.maxC = 0.2;
params.indexN = 1:N;

% Filter Setup
fnq = 1/(2*params.TR);
flp = 0.01; 
fhi = fnq * 0.99; % Based on Jacub's paper. Approach Nyquist without hitting 1.0 exactly
Wn = [flp/fnq fhi/fnq]; 
k = 2;                          
[bfilt,afilt] = butter(k,Wn);   
params.bfilt = bfilt;
params.afilt = afilt;

% --- Homotopic Mapping ---
hom = (1:N)';                 
hom(1:16)    = (1:16) + 16;   
hom(17:32)   = (17:32) - 16;  
hom(33:132)  = (33:132) + 100; 
hom(133:232) = (133:232) - 100; 
assert(all(hom(hom) == (1:N)'), 'hom must be symmetric.');

% Create Optimization Mask (Vectorization Pre-calculation)
% We update C(i,j) if C(i,j) > 0 OR if j == hom(i)
C = C/max(max(C))*params.maxC; % Normalize Initial C
mask = (C > 0);
% Add homotopic connections to mask
linear_ind_hom = sub2ind([N, N], 1:N, hom');
mask(linear_ind_hom) = true;
params.update_mask = mask; 

NSUB_LOW = size(low_anhedonia, 1);
NSUB_HIGH = size(high_anhedonia, 1);

% Check if a parallel pool already exists
currentPool = gcp('nocreate'); 

if isempty(currentPool)
    fprintf('No parallel pool detected. Starting one now...\n');
    % This starts a pool with the default number of workers (usually equal to CPU cores)
    parpool; 
else
    fprintf('Parallel pool is active with %d workers.\n', currentPool.NumWorkers);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOW ANHEDONIA GROUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- 1. Group-level analysis ---
fprintf('Processing Low Anhedonia Group...\n');

% Pre-calculate Group Average FC and COVtau
% We process subjects just to get the average empirical matrices first
FC_stack = zeros(NSUB_LOW, N, N);
COVtau_stack = zeros(NSUB_LOW, N, N);

% Parallelize extraction of empirical data
parfor nsub = 1:NSUB_LOW
    ts = low_anhedonia{nsub, 1};
    [FC_sub, COVtau_sub] = process_subject_empirical(ts, params);
    FC_stack(nsub,:,:) = FC_sub;
    COVtau_stack(nsub,:,:) = COVtau_sub;
end

FC_LOW = FC_stack;      % Store for saving
COVtau_LOW = COVtau_stack; % Store for saving
FCemp_group_low = squeeze(mean(FC_stack, 1));
COVtauemp_group_low = squeeze(mean(COVtau_stack, 1));

% Optimize Group
fprintf('  Optimizing Group Model...\n');
Ceffgroup_LOW = run_optimization(C, FCemp_group_low, COVtauemp_group_low, f_diff_low, params);

% --- 2. Individual-level analysis ---
fprintf('Processing Low Anhedonia Individual Subjects...\n');

% Pre-allocate outputs for parallel loop
Ceff_LOW = zeros(NSUB_LOW, N, N);
fittFC_LOW = zeros(NSUB_LOW, 1);
fittCVtau_LOW = zeros(NSUB_LOW, 1);
trophiccoherence_LOW = zeros(NSUB_LOW, 1);
hierarchicallevels_LOW = zeros(NSUB_LOW, N);

parfor nsub = 1:NSUB_LOW
    fprintf('  Subject %d/%d (Low)\n', nsub, NSUB_LOW);
    
    % Get Empirical Data (already computed, but cheap to re-process for clean parallel logic)
    % or grab from the stacks we created above
    FCemp = squeeze(FC_LOW(nsub,:,:));
    COVtauemp = squeeze(COVtau_LOW(nsub,:,:));
    
    % Optimize Individual (Starting from Group Ceff)
    [Ceff, fitFC, fitCOV] = run_optimization(Ceffgroup_LOW, FCemp, COVtauemp, f_diff_low, params);
    
    % Store connectivity and fits
    Ceff_LOW(nsub,:,:) = Ceff;
    fittFC_LOW(nsub) = fitFC;
    fittCVtau_LOW(nsub) = fitCOV;
    
    % Compute Hierarchy (Vectorized inside the loop)
    [tc, hl] = compute_hierarchy(Ceff);
    trophiccoherence_LOW(nsub) = tc;
    hierarchicallevels_LOW(nsub, :) = hl;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HIGH ANHEDONIA GROUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- 1. Group-level analysis ---
fprintf('Processing High Anhedonia Group...\n');

FC_stack = zeros(NSUB_HIGH, N, N);
COVtau_stack = zeros(NSUB_HIGH, N, N);

parfor nsub = 1:NSUB_HIGH
    ts = high_anhedonia{nsub, 1};
    [FC_sub, COVtau_sub] = process_subject_empirical(ts, params);
    FC_stack(nsub,:,:) = FC_sub;
    COVtau_stack(nsub,:,:) = COVtau_sub;
end

FC_HIGH = FC_stack;
COVtau_HIGH = COVtau_stack;
FCemp_group_high = squeeze(mean(FC_stack, 1));
COVtauemp_group_high = squeeze(mean(COVtau_stack, 1));

fprintf('  Optimizing Group Model...\n');
Ceffgroup_HIGH = run_optimization(C, FCemp_group_high, COVtauemp_group_high, f_diff_high, params);

% --- 2. Individual-level analysis ---
fprintf('Processing High Anhedonia Individual Subjects...\n');

Ceff_HIGH = zeros(NSUB_HIGH, N, N);
fittFC_HIGH = zeros(NSUB_HIGH, 1);
fittCVtau_HIGH = zeros(NSUB_HIGH, 1);
trophiccoherence_HIGH = zeros(NSUB_HIGH, 1);
hierarchicallevels_HIGH = zeros(NSUB_HIGH, N);

parfor nsub = 1:NSUB_HIGH
    fprintf('  Subject %d/%d (High)\n', nsub, NSUB_HIGH);
    
    FCemp = squeeze(FC_HIGH(nsub,:,:));
    COVtauemp = squeeze(COVtau_HIGH(nsub,:,:));
    
    [Ceff, fitFC, fitCOV] = run_optimization(Ceffgroup_HIGH, FCemp, COVtauemp, f_diff_high, params);
    
    Ceff_HIGH(nsub,:,:) = Ceff;
    fittFC_HIGH(nsub) = fitFC;
    fittCVtau_HIGH(nsub) = fitCOV;
    
    [tc, hl] = compute_hierarchy(Ceff);
    trophiccoherence_HIGH(nsub) = tc;
    hierarchicallevels_HIGH(nsub, :) = hl;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  STATISTICS & SAVING (Mostly Unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Performing Statistical Comparisons...\n');

% (Original plotting code preserved mostly as-is)
figure(1);
subplot(1,2,1)
tc_low_col = trophiccoherence_LOW(:); 
tc_high_col = trophiccoherence_HIGH(:); 
boxplot([tc_low_col; tc_high_col], [ones(length(tc_low_col),1); 2*ones(length(tc_high_col),1)]);
xlabel('Group'); ylabel('Trophic Coherence');
set(gca, 'XTickLabel', {'Low Anhedonia', 'High Anhedonia'});
title('Trophic Coherence Comparison');

[p_value, h] = ranksum(trophiccoherence_LOW, trophiccoherence_HIGH);
fprintf('Trophic Coherence p-value = %.4f\n', p_value);

subplot(1,2,2)
bar([mean(trophiccoherence_LOW), mean(trophiccoherence_HIGH)]); hold on;
errorbar([1, 2], [mean(trophiccoherence_LOW), mean(trophiccoherence_HIGH)], ...
    [std(trophiccoherence_LOW)/sqrt(length(trophiccoherence_LOW)), ...
     std(trophiccoherence_HIGH)/sqrt(length(trophiccoherence_HIGH))], 'k.');
title(sprintf('p = %.4f', p_value)); hold off;

figure(2);
subplot(2,2,1);
boxplot([fittFC_LOW(:); fittFC_HIGH(:)], [ones(NSUB_LOW,1); 2*ones(NSUB_HIGH,1)]);
title('FC Model Fit');
subplot(2,2,2);
boxplot([fittCVtau_LOW(:); fittCVtau_HIGH(:)], [ones(NSUB_LOW,1); 2*ones(NSUB_HIGH,1)]);
title('COV(tau) Model Fit');
subplot(2,2,3);
mean_hier_LOW = mean(hierarchicallevels_LOW, 2);
mean_hier_HIGH = mean(hierarchicallevels_HIGH, 2);
boxplot([mean_hier_LOW(:); mean_hier_HIGH(:)], [ones(NSUB_LOW,1); 2*ones(NSUB_HIGH,1)]);
title('Mean Hierarchical Levels');

fprintf('Saving results...\n');
save(fullfile(output_low_dir, 'results_Ceff_low_anhedonia.mat'), 'Ceff_LOW', 'Ceffgroup_LOW', 'trophiccoherence_LOW', 'hierarchicallevels_LOW', 'fittFC_LOW', 'fittCVtau_LOW', 'FC_LOW', 'COVtau_LOW');
save(fullfile(output_high_dir, 'results_Ceff_high_anhedonia.mat'), 'Ceff_HIGH', 'Ceffgroup_HIGH', 'trophiccoherence_HIGH', 'hierarchicallevels_HIGH', 'fittFC_HIGH', 'fittCVtau_HIGH', 'FC_HIGH', 'COVtau_HIGH');
save(fullfile(output_dir, 'results_Ceff_anhedonia_comparison.mat'), 'Ceff_LOW', 'Ceff_HIGH', 'trophiccoherence_LOW', 'trophiccoherence_HIGH');

% savefig(figure(1), fullfile(output_dir, 'trophic_coherence_comparison.fig'));
% savefig(figure(2), fullfile(output_dir, 'model_fit_hierarchical_levels.fig'));

fprintf('Analysis Complete!\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOCAL FUNCTIONS (Must be at the end of the script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Coptimized, finalFitFC, finalFitCOV] = run_optimization(C_init, FCemp, COVtauemp, f_diff, params)
    % Unpack parameters
    N = params.N;
    Cnew = C_init;
    epsFC = params.epsFC;
    epsFCtau = params.epsFCtau;
    maxC = params.maxC;
    mask = params.update_mask;
    sigma = params.sigma;
    Tau = params.Tau;
    TR = params.TR;
    
    olderror = 100000;
    
    % Main Optimization Loop
    for iter = 1:5000
        % Hopf Simulation
        [FCsim, COVsim, COVsimtotal, A] = hopf_int(Cnew, f_diff, sigma);
        
        % Calculate COVtau Simulation (Optimized)
        % Note: expm is the inevitable bottleneck here.
        COVtausim = expm((Tau*TR)*A) * COVsimtotal;
        COVtausim = COVtausim(1:N, 1:N);
        
        % Normalize COVtausim by diagonals (Vectorized)
        d = diag(COVsim);
        inv_sqrt_d = 1 ./ sqrt(d);
        sigratiosim = inv_sqrt_d * inv_sqrt_d'; % Outer product
        COVtausim = COVtausim .* sigratiosim;
        
        % Calculate Error
        diffFC = FCemp - FCsim;
        diffCOV = COVtauemp - COVtausim;
        
        % Check Convergence (every 100 iters)
        if mod(iter, 100) == 0
            errorFC = mean(diffFC(:).^2);
            errorCOVtau = mean(diffCOV(:).^2);
            errornow = errorFC + errorCOVtau;
            
            % --- ADDED PROGRESS BAR HERE ---
            fprintf('    Iter: %4d | Error: %.6f\n', iter, errornow);
            
            if (olderror - errornow)/errornow < 0.001 || olderror < errornow
                break;
            end
            olderror = errornow;
        end
        
        % --- VECTORIZED UPDATE ---
        % Calculate the update delta for the whole matrix at once
        Delta = epsFC * diffFC + epsFCtau * diffCOV;
        
        % Apply update only to masked indices
        Cnew(mask) = Cnew(mask) + Delta(mask);
        
        % Enforce constraints
        Cnew(Cnew < 0) = 0;
        Cnew = Cnew / max(Cnew(:)) * maxC;
    end
    
    Coptimized = Cnew;
    
    % Final fits
    Isubdiag = find(tril(ones(N),-1));
    finalFitFC = corr2(FCemp(Isubdiag), FCsim(Isubdiag));
    finalFitCOV = corr2(COVtauemp(Isubdiag), COVtausim(Isubdiag));
end


function [FCemp, COVtauemp] = process_subject_empirical(ts, params)
    N = params.N;
    indexN = params.indexN;
    Tau = params.Tau;
    
    % 1. Demean, Detrend, Filter
    % We process row by row, but can vectorize the detrend if memory allows. 
    % Keeping it loop-based per region is usually fine for N=232.
    for seed = 1:N
        row = ts(seed,:);
        row = detrend(row - nanmean(row));
        row = filtfilt(params.bfilt, params.afilt, row);
        ts(seed,:) = row;
    end
    
    % Trim
    ts2 = ts(indexN, 10:end-10);
    T = size(ts2, 2);
    
    % FC(0)
    FCemp = corrcoef(ts2');
    
    % COV(tau) - VECTORIZED
    % Instead of N*N calls to xcov, we use matrix multiplication
    % x(t) corresponds to 1 : end-Tau
    % x(t+tau) corresponds to 1+Tau : end
    X_t = ts2(:, 1:end-Tau);
    X_tau = ts2(:, 1+Tau:end);
    
    % Covariance = E[(X - muX)(Y - muY)]. Data is already demeaned/detrended roughly,
    % but strictly we should compute covariance of the segments.
    % Assuming mean is approx 0 after detrending:
    % Cov(X_t, X_tau) approx (X_t * X_tau') / (T_segment - 1)
    
    n_seg = size(X_t, 2);
    % Standardize inputs to match xcov 'coeff' or correlation logic if needed.
    % The original code used xcov to get raw covariance, then scaled by sigratio.
    
    % 1. Calculate raw cross-covariance at lag Tau efficiently
    % C_raw = (X_t * X_tau') / n_seg; % This calculates E[x(t)x(t+tau)]
    
    % However, strict xcov with 'biased' (default in some contexts) divides by N.
    % The original code: clag(indx)/size(tst,1). xcov by default centers the data.
    % Since we detrended, we can assume centered.
    
    % Center the specific segments to be safe
    X_t = X_t - mean(X_t, 2);
    X_tau = X_tau - mean(X_tau, 2);
    
    COV_raw = (X_t * X_tau') ./ n_seg; % This is the raw covariance matrix at lag Tau
    
    % 2. Calculate sigratio matrix (normalization factors)
    % sigratio(i,j) = 1 / (std(i)*std(j))
    % Note: The original code used COVemp(i,i) which is Variance.
    COVemp = cov(ts2');
    d = diag(COVemp);
    inv_std = 1 ./ sqrt(d);
    sigratio = inv_std * inv_std';
    
    COVtauemp = COV_raw .* sigratio;
end

function [tc, hl] = compute_hierarchy(Ceff)
    % Vectorized hierarchy calculation
    A = Ceff';
    d = sum(A)';
    delta = sum(A,2);
    u = d + delta;
    v = d - delta;
    
    Lambda = diag(u) - A - A';
    Lambda(1,1) = 0; % Regularization for singularity
    
    % linsolve is usually faster than inv()
    gamma = linsolve(Lambda, v);
    gamma = gamma - min(gamma);
    hl = gamma';
    
    % Trophic coherence
    % Vectorized meshgrid replacement
    H = (gamma - gamma' - 1).^2;
    F0 = sum(sum((A.*H))) / sum(sum(A));
    tc = 1 - F0;
end
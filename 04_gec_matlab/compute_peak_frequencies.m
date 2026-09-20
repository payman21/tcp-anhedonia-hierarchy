% this is the modified code to perform the torphic levels analysis using
% clusters that were created using only SHAPS (no factor analysis)
clear all;

% Paths
data_dir = '/Users/proghani/Documents/personal/code_experiments/neuro/phd_stuff/tcp_parcellations/data/anhedonia/NEW_4_factor_clustering_ipnybV4';
output_dir = '/Users/proghani/Documents/personal/local_thesis_proj/TCP_thesis_project/my_analysis/trophic_coherence_analysis/outputs/NEW_4_factor_clustering_ipnybV4';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

N=232;
indexN=1:N;

% Parameters of the data
TR=0.8; % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR); % Nyquist frequency
flp = 0.01; % lowpass frequency of filter (Hz), based on Jacub's paper
fhi = 0.6; % highpass (based on the TR, close to Nyquist frequency - Jakub's paper)
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
k=2; % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn); % construct the filter

% Process both high and low anhedonia groups
groups = {'high_anhedonia', 'low_anhedonia'};

for g = 1:length(groups)
    group_name = groups{g};
    fprintf('Processing %s...\n', group_name);

    % Load data
    loaded = load(fullfile(data_dir, [group_name '.mat']));
    group_data = loaded.(group_name);

    NSUB = size(group_data, 1);
    clear f_diff_sub FCemp;

    for sub=1:NSUB
        sub
        clear signal_filt Power_Areas;

        ts = group_data{sub,1};
        ts = ts(indexN,:);

        for seed=1:N
            ts(seed,:) = detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) = filtfilt(bfilt,afilt,ts(seed,:));
        end

        signal_filt = signal_filt(:,10:end-10);
        [Ns, Tmaxred] = size(signal_filt);
        TT = Tmaxred;
        Ts = TT*TR;
        freq = (0:TT/2-1)/Ts;
        nfreqs = length(freq);

        for seed=1:N
            pw = abs(fft(zscore(signal_filt(seed,:))));
            PowSpect = pw(1:floor(TT/2)).^2/(TT/TR);
            Power_Areas = gaussfilt(freq,PowSpect,0.005);
            [~,index] = max(Power_Areas);
            index = squeeze(index);
            f_diff_sub(sub,seed) = freq(index);
        end

        ts = zscore(ts,[],2);
        FCemp(sub,:,:) = corrcoef(ts');
    end

    f_diff = mean(f_diff_sub,1);
    FCemp = squeeze(nanmean(FCemp));

    save(fullfile(output_dir, ['empirical_' group_name '.mat']), 'FCemp', 'f_diff');
    fprintf('Saved empirical_%s.mat to: %s\n', group_name, output_dir);
end

fprintf('All analysis complete.\n');

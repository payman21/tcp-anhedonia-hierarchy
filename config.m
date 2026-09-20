function p = config()
%CONFIG  Filesystem locations for the MATLAB stages of the pipeline.
%
%   None of these paths are distributed with this repository (see README).
%   Edit the defaults below, or set the matching environment variables.

p.TCP_DATA = getenv_or('TCP_DATA', '/path/to/tcp_parcellations/data');
p.GEC      = getenv_or('TCP_GEC',  '/path/to/trophic_coherence_analysis');

% Cluster-assigned parcellated time series and the group structural connectome
p.CLUSTERS = fullfile(p.TCP_DATA, 'anhedonia', 'NEW_4_factor_clustering_ipnybV4');
p.SC_FILE  = fullfile(p.TCP_DATA, 'sc_matrix', ...
    'SC_schaefer200_tian_S2_7Networks_32fold_groupconnectome_3T_MNI152NLin2009cAsym_2mm.mat');

% Where stage-04 results are written
p.OUTPUT   = fullfile(p.GEC, 'outputs', 'NEW_4_factor_clustering_ipnybV4');

% Third-party MATLAB toolboxes required by the rendering scripts (stage 06).
% Install separately; see README.
p.RENDER_UTILS = getenv_or('TCP_RENDER_UTILS', '/path/to/render_utils');
end

function v = getenv_or(name, default)
v = getenv(name);
if isempty(v), v = default; end
end

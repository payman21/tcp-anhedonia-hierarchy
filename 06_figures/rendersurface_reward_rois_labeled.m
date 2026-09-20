function rendersurface_reward_rois_labeled(roi_values, rangemin, rangemax, inv, clmap, surfacetype)
% Renders mean hierarchical level values for 4 cortical reward ROIs, with ROI labels.
% Adapted from MODIFIED_rendersurface_dbs80.m (ML Kringelbach, June 2020)
%
% roi_values : 4-element vector [mPFC, OFC, ACC, dlPFC]
%   mPFC  - Default_PFC parcels (pCun excluded)
%   OFC   - Limbic_OFC + Cont_OFC parcels
%   ACC   - SalVentAttn_Med + Cont_Cing parcels
%   dlPFC - Cont_PFCl parcels
%
% rangemin, rangemax : limits for colorscheme (default: min/max of roi_values)
%
% inv:
%  0 colormap, interp (default)
%  1 flip colormap, interp
%  2 colormap, only three colours
%
% clmap : from othercolor, default 'RdBu7'
%
% surfacetype:
%  1 midthickness
%  2 inflated (default)
%  3 very inflated
%
% Example usage:
%   rendersurface_reward_rois_labeled([0.21, 0.17, 0.19, 0.22])
% to run:
%   data = load('path/to/your_4_roi_values.txt');
%   rendersurface_reward_rois_labeled(data)

% Paths relative to this script's location
scriptdir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptdir, 'render_utils', 'gifti-main'));
addpath(fullfile(scriptdir, 'render_utils', 'xmltree-main'));
addpath(fullfile(scriptdir, 'render_utils', 'subtightplot'));
addpath(fullfile(scriptdir, 'render_utils', 'othercolor'));

if ~exist('rangemin', 'var') || isempty(rangemin)
    rangemin = min(roi_values);
end

if ~exist('rangemax', 'var') || isempty(rangemax)
    rangemax = max(roi_values);
end

if ~exist('inv', 'var')
    inv = 0;
end

if ~exist('clmap', 'var')
    clmap = 'RdBu7';
end

if ~exist('surfacetype', 'var')
    surfacetype = 2; % default is inflated
end

disp(rangemin)
disp(rangemax)

% make space tight
make_it_tight = true;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight, clear subplot; end

% Surface meshes
basedir = fullfile(scriptdir, 'render_utils');
glassers_L   = gifti(fullfile(basedir, 'Glasser360.L.mid.32k_fs_LR.surf.gii'));
glassersi_L  = gifti(fullfile(basedir, 'Glasser360.L.inflated.32k_fs_LR.surf.gii'));
glassersvi_L = gifti(fullfile(basedir, 'Glasser360.L.very_inflated.32k_fs_LR.surf.gii'));
glassersf_L  = gifti(fullfile(basedir, 'Glasser360.L.flat.32k_fs_LR.surf.gii'));
glassers_R   = gifti(fullfile(basedir, 'Glasser360.R.mid.32k_fs_LR.surf.gii'));
glassersi_R  = gifti(fullfile(basedir, 'Glasser360.R.inflated.32k_fs_LR.surf.gii'));
glassersvi_R = gifti(fullfile(basedir, 'Glasser360.R.very_inflated.32k_fs_LR.surf.gii'));
glassersf_R  = gifti(fullfile(basedir, 'Glasser360.R.flat.32k_fs_LR.surf.gii'));

switch surfacetype
    case 1
        display_surf_left  = glassers_L;
        display_surf_right = glassers_R;
    case 2
        display_surf_left  = glassersi_L;
        display_surf_right = glassersi_R;
    case 3
        display_surf_left  = glassersvi_L;
        display_surf_right = glassersvi_R;
end

sl = display_surf_left;
sr = display_surf_right;

% Atlas label GIFTIs (Schaefer200 + Tian S2, cortical surface extracted)
atlasdir = fullfile(scriptdir, 'my_schaefer_tian_files');
label_L = gifti(fullfile(atlasdir, 'Schaefer200_Tian_S2.L.32k_fs_LR.label.gii'));
label_R = gifti(fullfile(atlasdir, 'Schaefer200_Tian_S2.R.32k_fs_LR.label.gii'));

% Functional GIFTI templates (reuse DBS80 files — same 32k vertex count)
vl = gifti(fullfile(basedir, 'dbs80_left.func.gii'));
vr = gifti(fullfile(basedir, 'dbs80_right.func.gii'));

% Initialise all vertices to 0 (grey background / unlabelled cortex)
vl.cdata(:) = 0;
vr.cdata(:) = 0;

% ROI label ID arrays derived from Schaefer200_Tian_S2 LabelTable
% Order: {mPFC, OFC, ACC, dlPFC}
%
% mPFC  : Default_PFC parcels, pCun excluded
% OFC   : Limbic_OFC + Cont_OFC parcels
% ACC   : SalVentAttn_Med + Cont_Cing parcels
% dlPFC : Cont_PFCl parcels
roi_label_ids_L = { ...
    [115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127], ...  % mPFC
    [87, 88, 97], ...                                                         % OFC
    [84, 85, 86, 104, 105], ...                                               % ACC
    [98, 99, 100, 101, 102] ...                                               % dlPFC
};

roi_label_ids_R = { ...
    [222, 223, 224, 225, 226, 227, 228, 229], ...                            % mPFC
    [191, 192, 193], ...                                                      % OFC
    [188, 189, 190, 210, 211], ...                                            % ACC
    [202, 203, 204, 205, 206, 207, 208] ...                                   % dlPFC
};

% Paint each ROI with its scalar value and record vertex indices for labels
roi_names      = {'mPFC', 'OFC', 'ACC', 'dlPFC'};
roi_vertices_L = cell(4,1);
roi_vertices_R = cell(4,1);

for r = 1:4
    mask_L = false(size(vl.cdata));
    mask_R = false(size(vr.cdata));
    for k = roi_label_ids_L{r}
        mask_L = mask_L | (label_L.cdata == k);
        vl.cdata(label_L.cdata == k) = roi_values(r);
    end
    for k = roi_label_ids_R{r}
        mask_R = mask_R | (label_R.cdata == k);
        vr.cdata(label_R.cdata == k) = roi_values(r);
    end
    roi_vertices_L{r} = find(mask_L);
    roi_vertices_R{r} = find(mask_R);
end


%% rendering

    % create figure
    hfig = figure;
    set(gcf, 'Position', [100, 100, 500, 500]);

    subplot(3,2,1); % left hemisphere lateral view
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(1) = patch(ax2, 'Faces', sl.faces, 'vertices', sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(-90, 0);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sl.vertices, sl.faces, vl.cdata);
    label_rois(ax2, sl.vertices, roi_vertices_L, roi_names);

    subplot(3,2,3); % left hemisphere medial view
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(1) = patch(ax2, 'Faces', sl.faces, 'vertices', sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(90, 0);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sl.vertices, sl.faces, vl.cdata);
    label_rois(ax2, sl.vertices, roi_vertices_L, roi_names);

    subplot(3,2,4); % right hemisphere lateral view
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(2) = patch(ax2, 'Faces', sr.faces, 'vertices', sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(-90, 0);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sr.vertices, sr.faces, vr.cdata);
    label_rois(ax2, sr.vertices, roi_vertices_R, roi_names);

    subplot(3,2,2); % right hemisphere medial view
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(2) = patch(ax2, 'Faces', sr.faces, 'vertices', sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(90, 0);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sr.vertices, sr.faces, vr.cdata);
    label_rois(ax2, sr.vertices, roi_vertices_R, roi_names);

    % flat maps
    sl = glassersf_L;
    sr = glassersf_R;

    subplot(3,2,5); % left hemisphere flat
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(2) = patch(ax2, 'Faces', sl.faces, 'vertices', sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(0, 90);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sl.vertices, sl.faces, vl.cdata);
    label_rois(ax2, sl.vertices, roi_vertices_L, roi_names);

    subplot(3,2,6); % right hemisphere flat
    ax2 = gca;
    axis(ax2, 'equal');
    axis(ax2, 'off');
    s(2) = patch(ax2, 'Faces', sr.faces, 'vertices', sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor', 'interp', 'EdgeColor', 'none');
    set(ax2, 'CLim', [rangemin rangemax]);
    view(0, 90);
    camlight;
    lighting gouraud;
    material dull;
    draw_roi_borders(ax2, sr.vertices, sr.faces, vr.cdata);
    label_rois(ax2, sr.vertices, roi_vertices_R, roi_names);
    colorbar('southoutside');

    switch inv
        case 0
            c = othercolor(clmap);
            colormap(c)
        case 1
            c = flipud(othercolor(clmap));
        case 2
            c = othercolor(clmap, 3);
    end
    % neutral grey for unlabelled cortex (background vertices stay at 0)
    c = othercolor(clmap);
    c(1,1) = 0.95; c(1,2) = 0.95; c(1,3) = 0.95;
    colormap(c)
    colorbar('Location', 'southoutside', 'Position', [.25 .05 .5 .05]);

end


function draw_roi_borders(ax, vertices, faces, vertex_data)
% Draw thin black outlines around the boundary of ROI regions.
%
% A "boundary edge" is an edge shared by exactly one ROI face — i.e. it
% sits on the border between an ROI region and the grey background.
% Interior edges shared by two ROI faces are not drawn.
%
% Uses a single line() call (NaN-separated segments) for speed.

face_in_roi = any(vertex_data(faces) ~= 0, 2);
roi_faces   = faces(face_in_roi, :);

if isempty(roi_faces)
    return;
end

% All edges of ROI faces; sort each edge so the smaller vertex index is first
edges = [roi_faces(:,[1,2]); roi_faces(:,[2,3]); roi_faces(:,[1,3])];
edges = sort(edges, 2);

% Boundary edges appear exactly once (interior edges shared by 2 ROI faces appear twice)
[uniq_edges, ~, ic] = unique(edges, 'rows');
counts   = accumarray(ic, 1);
boundary = uniq_edges(counts == 1, :);

if isempty(boundary)
    return;
end

% Build NaN-separated coordinate arrays — one line() call is much faster than a loop
n  = size(boundary, 1);
vx = [vertices(boundary(:,1),1), vertices(boundary(:,2),1), nan(n,1)]';
vy = [vertices(boundary(:,1),2), vertices(boundary(:,2),2), nan(n,1)]';
vz = [vertices(boundary(:,1),3), vertices(boundary(:,2),3), nan(n,1)]';

line(ax, vx(:), vy(:), vz(:), 'Color', 'k', 'LineWidth', 0.5);

end


function label_rois(ax, vertices, roi_vertex_idx, roi_names)
% Place a text label at the centroid of each ROI's vertex cloud.

for r = 1:length(roi_names)
    idx = roi_vertex_idx{r};
    if isempty(idx), continue; end
    cx = mean(vertices(idx, 1));
    cy = mean(vertices(idx, 2));
    cz = mean(vertices(idx, 3));
    text(ax, cx, cy, cz, roi_names{r}, ...
        'FontSize', 7, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Color', 'k', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'none', ...
        'Margin', 1);
end

end

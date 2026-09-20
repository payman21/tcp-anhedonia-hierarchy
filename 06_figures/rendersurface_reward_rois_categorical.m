function rendersurface_reward_rois_categorical(surfacetype)
% Renders 4 cortical reward ROIs each in a distinct categorical color.
% No scalar values needed — each ROI gets a fixed RGB color, making
% the 4 regions immediately distinguishable.
%
% ROIs and their colors:
%   mPFC  - blue   [Default_PFC parcels, pCun excluded]
%   OFC   - orange [Limbic_OFC + Cont_OFC parcels]
%   ACC   - green  [SalVentAttn_Med + Cont_Cing parcels]
%   dlPFC - red    [Cont_PFCl parcels]
%
% surfacetype:
%  1 midthickness
%  2 inflated (default)
%  3 very inflated
%
% Example usage:
%   rendersurface_reward_rois_categorical()
%   rendersurface_reward_rois_categorical(1)

scriptdir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptdir, 'render_utils', 'gifti-main'));
addpath(fullfile(scriptdir, 'render_utils', 'xmltree-main'));
addpath(fullfile(scriptdir, 'render_utils', 'subtightplot'));

if ~exist('surfacetype', 'var')
    surfacetype = 2;
end

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

% Atlas label GIFTIs
atlasdir = fullfile(scriptdir, 'my_schaefer_tian_files');
label_L = gifti(fullfile(atlasdir, 'Schaefer200_Tian_S2.L.32k_fs_LR.label.gii'));
label_R = gifti(fullfile(atlasdir, 'Schaefer200_Tian_S2.R.32k_fs_LR.label.gii'));

% Functional GIFTI templates (reuse DBS80 files — same 32k vertex count)
vl = gifti(fullfile(basedir, 'dbs80_left.func.gii'));
vr = gifti(fullfile(basedir, 'dbs80_right.func.gii'));

nvert_L = numel(vl.cdata);
nvert_R = numel(vr.cdata);

% Colorblind-friendly categorical colors: mPFC, OFC, ACC, dlPFC
roi_colors = [
    0.267, 0.467, 0.667;   % mPFC  — blue
    0.933, 0.467, 0.200;   % OFC   — orange
    0.133, 0.533, 0.200;   % ACC   — green
    0.800, 0.200, 0.067;   % dlPFC — red
];
roi_names = {'mPFC', 'OFC', 'ACC', 'dlPFC'};

% Background colour for unlabelled cortex
grey = [0.85, 0.85, 0.85];

% Initialise all vertices to grey (Nx3 RGB arrays)
rgb_L = repmat(grey, nvert_L, 1);
rgb_R = repmat(grey, nvert_R, 1);

% Scalar ROI index per vertex (0 = background) — used only for border detection
roi_idx_L = zeros(nvert_L, 1);
roi_idx_R = zeros(nvert_R, 1);

% ROI label ID arrays derived from Schaefer200_Tian_S2 LabelTable
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

roi_vertices_L = cell(4, 1);
roi_vertices_R = cell(4, 1);

for r = 1:4
    mask_L = false(nvert_L, 1);
    mask_R = false(nvert_R, 1);
    for k = roi_label_ids_L{r}
        mask_L = mask_L | (label_L.cdata == k);
    end
    for k = roi_label_ids_R{r}
        mask_R = mask_R | (label_R.cdata == k);
    end
    rgb_L(mask_L, :) = repmat(roi_colors(r, :), sum(mask_L), 1);
    rgb_R(mask_R, :) = repmat(roi_colors(r, :), sum(mask_R), 1);
    roi_idx_L(mask_L) = r;
    roi_idx_R(mask_R) = r;
    roi_vertices_L{r} = find(mask_L);
    roi_vertices_R{r} = find(mask_R);
end


%% Rendering

figure;
set(gcf, 'Position', [100, 100, 500, 500]);

subplot(3,2,1); % left hemisphere lateral view
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sl.faces, 'Vertices', sl.vertices, ...
    'FaceVertexCData', rgb_L, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sl.vertices, sl.faces, roi_idx_L);
label_rois(ax, sl.vertices, roi_vertices_L, roi_names, roi_colors);

subplot(3,2,3); % left hemisphere medial view
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sl.faces, 'Vertices', sl.vertices, ...
    'FaceVertexCData', rgb_L, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(90, 0); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sl.vertices, sl.faces, roi_idx_L);
label_rois(ax, sl.vertices, roi_vertices_L, roi_names, roi_colors);

subplot(3,2,2); % right hemisphere medial view
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sr.faces, 'Vertices', sr.vertices, ...
    'FaceVertexCData', rgb_R, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(90, 0); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sr.vertices, sr.faces, roi_idx_R);
label_rois(ax, sr.vertices, roi_vertices_R, roi_names, roi_colors);

subplot(3,2,4); % right hemisphere lateral view
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sr.faces, 'Vertices', sr.vertices, ...
    'FaceVertexCData', rgb_R, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(-90, 0); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sr.vertices, sr.faces, roi_idx_R);
label_rois(ax, sr.vertices, roi_vertices_R, roi_names, roi_colors);

% Flat maps
sl_flat = glassersf_L;
sr_flat = glassersf_R;

subplot(3,2,5); % left hemisphere flat
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sl_flat.faces, 'Vertices', sl_flat.vertices, ...
    'FaceVertexCData', rgb_L, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(0, 90); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sl_flat.vertices, sl_flat.faces, roi_idx_L);
label_rois(ax, sl_flat.vertices, roi_vertices_L, roi_names, roi_colors);

subplot(3,2,6); % right hemisphere flat
ax = gca; axis(ax, 'equal'); axis(ax, 'off');
patch(ax, 'Faces', sr_flat.faces, 'Vertices', sr_flat.vertices, ...
    'FaceVertexCData', rgb_R, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(0, 90); camlight; lighting gouraud; material dull;
draw_roi_borders(ax, sr_flat.vertices, sr_flat.faces, roi_idx_R);
label_rois(ax, sr_flat.vertices, roi_vertices_R, roi_names, roi_colors);

% Legend — invisible NaN patches added to last axes
hold(ax, 'on');
hp = gobjects(4, 1);
for r = 1:4
    hp(r) = patch(ax, [NaN NaN NaN], [NaN NaN NaN], [NaN NaN NaN], ...
        'FaceColor', roi_colors(r, :), 'EdgeColor', 'none');
end
legend(ax, hp, roi_names, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off', 'FontSize', 10, 'FontWeight', 'bold');

end


function draw_roi_borders(ax, vertices, faces, roi_idx)
% Draw black outlines around ROI boundaries.
% roi_idx is a scalar per vertex: 0 = background, 1-4 = ROI.

face_in_roi = any(roi_idx(faces) ~= 0, 2);
roi_faces   = faces(face_in_roi, :);
if isempty(roi_faces), return; end

edges = [roi_faces(:,[1,2]); roi_faces(:,[2,3]); roi_faces(:,[1,3])];
edges = sort(edges, 2);
[uniq_edges, ~, ic] = unique(edges, 'rows');
counts   = accumarray(ic, 1);
boundary = uniq_edges(counts == 1, :);
if isempty(boundary), return; end

n  = size(boundary, 1);
vx = [vertices(boundary(:,1),1), vertices(boundary(:,2),1), nan(n,1)]';
vy = [vertices(boundary(:,1),2), vertices(boundary(:,2),2), nan(n,1)]';
vz = [vertices(boundary(:,1),3), vertices(boundary(:,2),3), nan(n,1)]';
line(ax, vx(:), vy(:), vz(:), 'Color', 'k', 'LineWidth', 0.5);

end


function label_rois(ax, vertices, roi_vertex_idx, roi_names, roi_colors)
% Place a text label at each ROI's centroid, coloured to match the ROI.

for r = 1:length(roi_names)
    idx = roi_vertex_idx{r};
    if isempty(idx), continue; end
    cx = mean(vertices(idx, 1));
    cy = mean(vertices(idx, 2));
    cz = mean(vertices(idx, 3));
    text(ax, cx, cy, cz, roi_names{r}, ...
        'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Color', roi_colors(r, :), ...
        'BackgroundColor', 'white', ...
        'EdgeColor', roi_colors(r, :), ...
        'Margin', 2);
end

end

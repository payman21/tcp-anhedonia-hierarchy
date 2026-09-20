function rendersurface_schaefer200(schaefer200vector, rangemin, rangemax, inv, clmap, surfacetype)
% Renders a Schaefer 200 cortical parcellation vector onto brain surfaces.
% Adapted from rendersurface_dbs80.m (ML Kringelbach, June 2020)
%
% schaefer200vector : 200-element vector of values to render
%   indices   1:100  -> left  hemisphere cortical parcels (label IDs 33-132)
%   indices 101:200  -> right hemisphere cortical parcels (label IDs 133-232)
%   (matches the ordering in Schaefer200_Tian_S2 label files)
%
% rangemin, rangemax : limits for colorscheme
%
% inv:
%  0 colormap, interp
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
%   data = load('path/to/cortical_hl_low_rendering.txt');
%   rendersurface_schaefer200(data)

% Paths relative to this script's location
scriptdir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptdir, 'render_utils', 'gifti-main'));
addpath(fullfile(scriptdir, 'render_utils', 'xmltree-main'));
addpath(fullfile(scriptdir, 'render_utils', 'subtightplot'));
addpath(fullfile(scriptdir, 'render_utils', 'othercolor'));

if ~exist('rangemin', 'var')
    rangemin = min(schaefer200vector);
end

if ~exist('rangemax', 'var')
    rangemax = max(schaefer200vector);
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

% Initialise all vertices to 0 (background / unlabelled medial wall)
vl.cdata(:) = 0;
vr.cdata(:) = 0;

% Map left hemisphere values:
% schaefer200vector(1:100)   -> label IDs 33:132
for i = 1:100
    label_id = 32 + i;
    idx = find(label_L.cdata == label_id);
    vl.cdata(idx) = schaefer200vector(i);
end

% Map right hemisphere values:
% schaefer200vector(101:200) -> label IDs 133:232
for i = 1:100
    label_id = 132 + i;
    idx = find(label_R.cdata == label_id);
    vr.cdata(idx) = schaefer200vector(100 + i);
end


%% rendering

    % create figure
    hfig = figure;
    set(gcf, 'Position',  [100, 100, 500, 500]);
    
    subplot(3,2,1); %left hemisphere side view
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0);
    camlight;
    lighting gouraud;
    material dull;


    subplot(3,2,3); %left hemisphere midline
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    subplot(3,2,4); %right hemisphere side view
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0)
    camlight;
    lighting gouraud;
    material dull;
 
    subplot(3,2,2); %right hemisphere midline
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    % flatmaps (commented out — only showing the 4 inflated views)
    % sl=glassersf_L;
    % sr=glassersf_R;
    %
    % subplot(3,2,5); %left hemisphere flat
    % ax2=gca;
    % axis(ax2,'equal');
    % axis(ax2,'off');
    % s(2) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    % set(ax2,'CLim',[rangemin rangemax]);
    % view(0,90)
    % camlight;
    % lighting gouraud;
    % material dull;
    % %colorbar('southoutside');
    %
    % subplot(3,2,6); %right hemisphere flat
    % ax2=gca;
    % axis(ax2,'equal');
    % axis(ax2,'off');
    % s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    % set(ax2,'CLim',[rangemin rangemax]);
    % view(0,90)
    % camlight;
    % lighting gouraud;
    % material dull;
    % colorbar('southoutside');
    
    switch inv
        case 0
            % use selected colormap (interpolated to 64 values)
            c=othercolor(clmap);
            colormap(c)  
            disp(c)
        case 1
            % flip colormap (interpolated to 64 values)
            c=flipud(othercolor(clmap));
        case 2
            % use specialised version with only three values
            c=othercolor(clmap,3);
    end;
    % neutral brain colour: grey
    c=othercolor(clmap);
    c(1,1)=0.95;c(1,2)=0.95;c(1,3)=0.95; 
    colormap(c)
    colorbar('Location', 'southoutside', 'Position', [.25 .34 .5 .05]);

end


%% render_brainnet.m
% Schematic glass-brain rendering of prefrontal + striatal ROIs.
%
% Renders the labelled volume built by make_rois.py as coloured 3D blobs
% inside a semi-transparent ICBM152 surface, one PNG per viewing angle.
%
% Run make_rois.py first. Then just press Run -- no arguments, no GUI
% interaction. Output lands in output/.
%
% Everything you would normally want to change is in the CONFIG block.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
if isempty(here)          % pasted into the command window rather than Run
    here = pwd;
end

bnv = fullfile(here, 'BrainNet-Viewer');
addpath(bnv);
addpath(fullfile(bnv, 'BrainNet_spm_files'));   % supplies spm_vol/spm_read_vols

%% ------------------------------- CONFIG --------------------------------

% Which ROI sets to render. Each is built by make_rois.py:
%   'split'    -> 9 regions, caudate/putamen/NAcc drawn separately
%   'striatum' -> 7 regions, those three merged into one blob
%   'striatal' -> 3 regions, striatum only, no cortex
% Colours are consistent across all three -- caudate is the same orange
% everywhere.
variants = {'split', 'striatal'};

meshAlpha = 0.13;                 % surface transparency; lower = glassier
meshColor = [0.92 0.92 0.92];     % surface tint
bakColor  = [1 1 1];              % background

% BrainMesh_ICBM152_smoothed is the cleaner-looking of the two ICBM meshes.
% Drop the _smoothed for a more gyrified surface.
surfFile = fullfile(bnv, 'Data', 'SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');

imgWidth  = 2400;
imgHeight = 1800;
imgDPI    = 600;

% name, azimuth, elevation. The brain is transparent, so a lateral view
% already shows the medial ROIs (sgACC, pACC, vmPFC) through the near
% hemisphere -- there is no separate "medial" view to render.
views = { ...
    'left_lateral',  -90,  0 ; ...
    'right_lateral',  90,  0 ; ...
    'anterior',        0,  0 ; ...
    'dorsal',          0, 90 };

renderContactSheet = true;   % extra 8-panel overview image

% Bare glass brain with no ROIs at all, for use as an underlay or a
% "regions omitted" panel. Names must match the first column of `views`;
% {} to skip. Same camera and canvas as the ROI renders, so the pair
% overlays pixel-for-pixel.
emptyViews = {'anterior', 'left_lateral'};

%% ------------------------------------------------------------------------

outDir = fullfile(here, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% Start from BrainNet's own defaults so the saved cfg is complete.
% BrainNet_MapCfg only backfills *some* missing EC fields, so hand-rolling
% the struct from scratch is fragile. Opening the GUI once populates the
% global EC with every default; we snapshot it and close.
h0 = BrainNet;
global EC
base = EC;
close(h0);

%% Settings shared by every render
base.msh.alpha      = meshAlpha;
base.msh.color      = meshColor;
base.msh.color_type = 1;          % single colour, not per-vertex
base.bak.color      = bakColor;

base.glb.material       = 'dull';
base.glb.shading        = 'interp';
base.glb.lighting       = 'gouraud';   % 'phong' is deprecated in modern MATLAB
base.glb.lightdirection = 'right';
base.glb.detail         = 3;
base.glb.lr             = 0;           % 1 draws L/R letters onto the render

base.img.width  = imgWidth;
base.img.height = imgHeight;
base.img.dpi    = imgDPI;

% view_direction must be 4, otherwise BrainNet ignores view_az/view_el and
% falls back to a hard-coded left lateral.
base.lot.view = 1;
base.lot.view_direction = 4;

%% Render each ROI set
for iv = 1:numel(variants)
    variant = variants{iv};

    volFile = fullfile(here, 'rois', sprintf('schematic_rois_%s.nii', variant));
    colFile = fullfile(here, 'rois', sprintf('schematic_rois_%s_colors.txt', variant));
    cfgFile = fullfile(outDir, sprintf('cfg_%s.mat', variant));

    assert(exist(volFile, 'file') == 2, 'Missing %s -- run make_rois.py first.', volFile);
    assert(exist(colFile, 'file') == 2, 'Missing %s -- run make_rois.py first.', colFile);

    colors = load(colFile);
    nROI   = size(colors, 1);
    fprintf('%s: %d regions\n', variant, nROI);

    cfg = base;

    % Volume -> ROI drawing
    cfg.vol.type    = 2;         % 2 = draw ROIs, 1 = map values onto surface
    cfg.vol.display = 1;
    cfg.vol.pn = 0;  cfg.vol.nn = 0;
    cfg.vol.nx = 0;  cfg.vol.px = nROI;

    % drawall is only consulted by the Option GUI. In command-line mode the
    % list of labels to draw has to be given explicitly or nothing is drawn.
    cfg.vol.roi.drawall = 0;
    cfg.vol.roi.draw    = (1:nROI)';

    % Colour row i applies to the i-th entry of roi.draw, so row order here
    % must match the label numbering make_rois.py assigned. It does -- both
    % come from the same REGIONS list.
    cfg.vol.roi.color(1:nROI, :) = colors;
    cfg.vol.roi.colort           = cfg.vol.roi.color;

    cfg.vol.roi.smooth        = 1;
    cfg.vol.roi.smooth_kernal = 1;   % 1 = box, 2 = Gaussian
    cfg.vol.threshold   = 0;
    cfg.vol.clustersize = 0;

    for v = 1:size(views, 1)
        cfg.lot.view = 1;
        cfg.lot.view_az = views{v, 2};
        cfg.lot.view_el = views{v, 3};

        S.EC = cfg;
        save(cfgFile, '-struct', 'S');

        outPng = fullfile(outDir, sprintf('%s_%s.png', variant, views{v, 1}));
        fprintf('  rendering %-15s -> %s\n', views{v, 1}, outPng);

        h = BrainNet_MapCfg(surfFile, volFile, cfgFile, outPng);
        close(h);
    end

    if renderContactSheet
        cfg.lot.view = 2;
        S.EC = cfg;
        save(cfgFile, '-struct', 'S');

        outPng = fullfile(outDir, sprintf('%s_contactsheet.png', variant));
        fprintf('  rendering %-15s -> %s\n', 'contact sheet', outPng);

        h = BrainNet_MapCfg(surfFile, volFile, cfgFile, outPng);
        close(h);
    end
end

%% Bare glass brain, no ROIs
% Passing no volume file leaves FLAG.Loadfile == 1, so BrainNet draws the
% mesh only and ignores everything under EC.vol.
cfg = base;
cfgFile = fullfile(outDir, 'cfg_emptybrain.mat');

for k = 1:numel(emptyViews)
    v = find(strcmp(views(:, 1), emptyViews{k}), 1);
    assert(~isempty(v), 'emptyViews entry "%s" is not in `views`.', emptyViews{k});

    cfg.lot.view_az = views{v, 2};
    cfg.lot.view_el = views{v, 3};

    S.EC = cfg;
    save(cfgFile, '-struct', 'S');

    outPng = fullfile(outDir, sprintf('emptybrain_%s.png', views{v, 1}));
    fprintf('  rendering %-15s -> %s\n', ['empty ' views{v, 1}], outPng);

    h = BrainNet_MapCfg(surfFile, cfgFile, outPng);
    close(h);
end

fprintf('\nDone. Images in %s\n', outDir);

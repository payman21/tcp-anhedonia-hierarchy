function rendersurface_dbs80(dbs80vector,rangemin,rangemax, inv,clmap,surfacetype)
% script for rendering a dbs80 vector
%   ML Kringelbach June 2020
%
% dk68vector : the dbs80 values to be rendered
%
% rangemin, rangemax : limits for colorscheme
%
% inv:
%  0 colormap, interp
%  1 flip colormap, interp
%  2 colormap, only three colours
%
% clmap : from othercolor, default 
%
% surfacetype: 
%  1 midthickness
%  2 inflated (default)
%  3 very inflated
%addpath('D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\gifti-main\')
addpath('/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/gifti-main/')
%addpath('D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\xmltree-main\')
addpath('/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/xmltree-main/')


if ~exist('rangemin','var')
     rangemin=min(dbs80vector);
end

if ~exist('rangemax','var')
     rangemax=max(dbs80vector);
end

if ~exist('inv','var')
     inv=max(dbs80vector);
end

if ~exist('clmap','var')
    clmap= 'RdBu7';%'RdBu7'; %GnBu7
end

if ~exist('surfacetype','var')
     surfacetype=2; % default is inflated
end
disp(rangemin)
disp(rangemax)
%addpath('D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\subtightplot')
addpath('/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/subtightplot')
%addpath('D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\othercolor')
addpath('/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/othercolor')
% make space tight
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

% load the different views
% base='/Users/mortenk/Documents/MATLAB/osl/std_masks/';
% display_surf_left=gifti([base 'ParcellationPilot.L.inflated.32k_fs_LR.surf.gii']);
% display_surf_right=gifti([base 'ParcellationPilot.R.inflated.32k_fs_LR.surf.gii']);

%basedir='D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\';
basedir='/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/';
glassers_L=gifti([basedir 'Glasser360.L.mid.32k_fs_LR.surf.gii']);
glassersi_L=gifti([basedir 'Glasser360.L.inflated.32k_fs_LR.surf.gii']);
glassersvi_L=gifti([basedir 'Glasser360.L.very_inflated.32k_fs_LR.surf.gii']);
glassersf_L=gifti([basedir 'Glasser360.L.flat.32k_fs_LR.surf.gii']);
glassers_R=gifti([basedir 'Glasser360.R.mid.32k_fs_LR.surf.gii']);
glassersi_R=gifti([basedir 'Glasser360.R.inflated.32k_fs_LR.surf.gii']);
glassersvi_R=gifti([basedir 'Glasser360.R.very_inflated.32k_fs_LR.surf.gii']);
glassersf_R=gifti([basedir 'Glasser360.R.flat.32k_fs_LR.surf.gii']);

switch surfacetype
    case 1
        display_surf_left=glassers_L;
        display_surf_right=glassers_R;
    case 2
        display_surf_left=glassersi_L;
        display_surf_right=glassersi_R;
    case 3
        display_surf_left=glassersvi_L;
        display_surf_right=glassersvi_R;
end;

sl = display_surf_left;
sr = display_surf_right;



%base='D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\';
base='/media/irene/D450-97BC/Irene/Trophic_soldiers/Scripts/renders/render_utils/';
% label gifti
label_L=gifti([base 'fsaverage.L.dbs80_Atlas.32k_fs_LR.label.gii']);
label_R=gifti([base 'fsaverage.R.dbs80_Atlas.32k_fs_LR.label.gii']);
% functional gifti
%base='D:\Irene\Trophic_soldiers\Scripts\renders\render_utils\';
vl=gifti([ base 'dbs80_left.func.gii']);
vr=gifti([ base 'dbs80_right.func.gii']);

mindboggle=[2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,34,35];
for i=1:31
    idx{i}=find(label_L.cdata==mindboggle(i));
    vl.cdata(idx{i})=dbs80vector(i);
end;

% replace right hemisphere labels with correct values from
% dbsvector(80:-1:50)
for i=1:31
    idx{i}=find(label_R.cdata==mindboggle(i));
    vr.cdata(idx{i})=dbs80vector(81-i);
end;

% remove -1 from labels
idx{i}=find(label_L.cdata==-1);
vl.cdata(idx{i})=0;
idx{i}=find(label_R.cdata==-1);
vr.cdata(idx{i})=0;


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

    % flatmaps
    sl=glassersf_L;
    sr=glassersf_R;
    
    subplot(3,2,5); %left hemisphere flat
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    %colorbar('southoutside');

    subplot(3,2,6); %right hemisphere flat
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    colorbar('southoutside');
    
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
    colorbar('Location', 'southoutside', 'Position', [.25 .05 .5 .05]); 

end

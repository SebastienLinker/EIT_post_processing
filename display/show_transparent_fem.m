function [ h_fem, h_tr ] = show_transparent_fem( mdl, varargin )
%SHOW_TRANSPARENT_FEM Similar to show_fem, but adds the possibility to add
%transparency in the FE model
%   Addition of transparency should help to plot both reconstruction and
%   target at the same time

%	(C) 2015/05/15 Sebastien Martin

[trans_th, opts] = process_args(varargin{:});

if strcmp(mdl.type,'inv_model'); mdl = mdl.fwd_model; end;
if strcmp(mdl.type,'fwd_model')
    h_fem = show_fem(mdl, opts); h_tr = []; return;
else
    h_fem = show_fem(mdl.fwd_model, opts); % Show forward model only
end
% 2D models, don't need transparency
if size(mdl.fwd_model.nodes,2)==2;
    eidors_msg('@@@ 2D model, transparency is available on 3D models only',4);
    h_fem = show_fem(mdl, opts); h_tr = []; return;
end

nodes = mdl.fwd_model.nodes;
elems = mdl.fwd_model.elems;
colours = get_img_data(mdl);
h_tr = repaint_inho_local( colours, 'use_global', nodes, elems, [], mdl );
arrayfun(@(x) set(x, 'FaceAlpha', trans_th), h_tr);

% Show electrodes after repainting
show_electrodes_3d( mdl.fwd_model );

end

function [trans_th, opts] = process_args(varargin)
if nargin<2; opts = [0 0 0]; else opts=varargin{2}; end
if nargin<1; trans_th = .2; else trans_th=varargin{1}; end
end

%% Show conductivity
% Because EIDORS function repaint_homo does not output any handle, here we
% copy paste it
function hh = repaint_inho_local( mat,mat_ref,vtx,simp, thresh, clr_def )
if nargin<5
    thresh = [];
end
if nargin<6
    clr_def = [];
end
if strcmp(mat_ref, 'use_global')
    img.calc_colours.ref_level = mat_ref;
end

if isempty(thresh)
    thresh = 1/4;
end

% This line correctly handles node data
% patch('Vertices',vtx,'Faces',simp,'FaceVertexCData',colours','FaceColor','interp');

% looks best if eidors_colours.greylev < 0
[colours,scl_data] = calc_colours( mat, clr_def, 0);
ii=find( abs(scl_data) > thresh);
this_x = simp(ii,:);

colours= permute(colours(ii,:,:),[2,1,3]);
ELEM= vtx';

Xs=   zeros(3,length(ii));
Ys=   zeros(3,length(ii));
Zs=   zeros(3,length(ii));
switch(size(this_x,2))
    case 3
        idx_ = [1;2;3];
    case 4
        idx_ = [[1;2;3], ...
            [1;2;4], ...
            [1;3;4], ...
            [2;3;4]];
end
hh = [];
for idx=idx_
    Xs(:)=vtx(this_x(:,idx)',1);
    Ys(:)=vtx(this_x(:,idx)',2);
    Zs(:)=vtx(this_x(:,idx)',3);
    
    if size(colours,1)==1 && size(colours,2)==3
        % need to work around ^%$#%$# matlab bug which
        % forces an incorrect interpretation is colours of this size
        hh= [hh patch(Xs(:,[1:3,1]), ...
            Ys(:,[1:3,1]), ...
            Zs(:,[1:3,1]), ...
            colours(:,[1:3,1]), ...
            'EdgeColor','none','CDataMapping','direct')];
    else
        hh= [hh patch(Xs,Ys,Zs,colours, ...
            'EdgeColor','none','CDataMapping','direct')];
    end
end
end

%% Show the electrodes
% Use function show_fem-->show_electrodes_3d to do that
function show_electrodes_3d(mdl)
% show electrode positions on model
if ~isfield(mdl,'electrode'); return; end

ee= get_boundary( mdl );
for e=1:length(mdl.electrode)
    colour= electr_colour( e);
    if isfield(mdl.electrode(e),'pos') && ~isfield(mdl.electrode(e),'nodes')
        error('This case is not implemented yet');
    end
    elec_nodes= mdl.electrode(e).nodes;
    
    if length(elec_nodes) == 1  % point electrode model
        vtx= mdl.nodes(elec_nodes,:);
        line(vtx(1),vtx(2),vtx(3), ...
            'Marker','o','MarkerSize',12, ...
            'MarkerFaceColor',colour, 'MarkerEdgeColor', colour);
    else
        % find elems on boundary attached to this electrode
        map = zeros(length(mdl.nodes),1);
        map(mdl.electrode(e).nodes) = 1;
        ec = map(ee);
        sels= find(all(ec'));
        
        ee= get_boundary( mdl );
        paint_electrodes(sels,ee,mdl.nodes,colour);
    end
end

end

function ee= get_boundary( mdl )
if isfield(mdl,'boundary')
    ee= mdl.boundary;
else
    % calc and cache boundary
    ee = find_boundary( mdl.elems );
end
end

function colour= electr_colour( e);
if e==1;
    colour = [0,.7,0]; % light green electrode #1
elseif e==2
    colour = [0,.5,0]; % mid-green electrode #2
else
    colour = [0,.3,0]; % dark green
end
end

function paint_electrodes(sel,srf,vtx, colour);
if isempty(sel); return; end  % Not required after matlab 2014

[u n m] = unique(srf(sel,:));
fv.vertices = vtx(u,:);
fv.faces = reshape(m,[],3);
h = patch(fv,'FaceColor',colour);
% need 'direct' otherwise colourmap is screwed up
set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );
end
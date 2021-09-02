function [ conds, info_obj ] = createMultiEllipses( n_objs, n_imgs, overlap, mdl, bk, varargin )
%CREATEMULTIELLIPSES Creates multiple ellipsoids shapes in one single FEM
%   Use this file to test a specific case in which there are more than one
%   object to detect
%	@input
%		overlap: Whether or not the targets may overlap
%		varargin: The rest of the parameters to generate the ellipses
%	@output conds: The list of conductivities of each element in each image
%		info_obj: array containing size, position(X,Y,Z), rotation angles, and conductivity of
%		random objects (having the highest and lowest conductivities):
%		Total size= n_imgs*12
%
%	Works for 2D applications. TODO: 3D


debug = false;

try	mdl2=mdl.fwd_model; catch; mdl2=mdl; end; mdl=mdl2;
ell_props = parse_varargin(n_objs, varargin{:}); % Ellipses properties

n_elems = size(mdl.elems,1);
is2D = (size(mdl.nodes)==2);

conds = zeros(n_elems,n_imgs);

for k=1:1:n_imgs
    c_cond = zeros(n_elems,1);
    l=1; prev_objs = [];
    while (l<=n_objs)
        [new_cond, new_obj{l}] = createEllipseGeneral(1,mdl,bk,ell_props{:,l});
        elems_obj = find( new_cond(:)~=bk ); % Get elements representing the object
        if ~overlap & sum( ismember(elems_obj,prev_objs)>0 ) % Check overlapping
            continue;
        else
            prev_objs = [prev_objs; elems_obj];
        end
        new_cond = new_cond-bk;
        % Here, it depends on how do you want to deal with overlapping
        % 		c_cond( elems_obj ) = c_cond(elems_obj)+new_cond(elems_obj);
        c_cond( elems_obj ) = new_cond(elems_obj);
        
        if(l>1)
            if (new_obj{l}(4)>info_obj(4,k))
                info_obj([1:6],k) = new_obj{l};
            end
            if (new_obj{l}(4)<info_obj(8,k))
                info_obj([7:12],k) = new_obj{l};
            end
        else
            info_obj(:,k) = [new_obj{1} new_obj{1}];
        end
        if debug
            figure; show_fem( mk_image(mdl, c_cond) );
        end
        l=l+1;
    end
    conds(:,k) = c_cond + bk;
    if debug
        figure; show_fem( mk_image(mdl, conds(:,k)), 1 );
    end
end

end

% Different properties for each object
% Input argument for createEllipseGeneral function
% max_sz: maximal size of the object
% 			can be a 2*1 or 3*1 array with different sizes for 2 or 3 radii, default: 0.4
% min_sz: minimum size
%			can be a 2*1 or 3*1 array with different sizes for 2 or 3 radii, default: 0
% rad_probe: radius of the probe
% max_pos: maximal distance from the centre (default 1)
%				3D unbounded: We may need s pecific Z-axis: In this case,
%				please use a cell array (2*1) or a 3*1 vector with:
%				Maximal distance (default 1)
%				Position of Z-axis (minimum and maximum)
%				For fixed Z-axis, a 2*1 vector is also ok
% rot_angle: The rotation angle (2*1 array with minimum and maximum, in degrees, default [0 180])
% pos_angle: Position angle, (2*1 array with minimum and maximum, in degrees, default [0 360])
% cond_r
function [argout] = parse_varargin( n_objs, varargin )
% default values
max_sz=0.4; min_sz=0; rad_probe=0; max_pos=1; rot_angle=[0 180];
pos_angle=[0 360]; cond_r=[0 1];
switch length(varargin)
    case 1
        max_sz = varargin{1};
    case 2
        max_sz = varargin{1}; min_sz = varargin{2};
    case 3
        max_sz = varargin{1}; min_sz = varargin{2}; rad_probe = varargin{3};
    case 4
        max_sz = varargin{1}; min_sz = varargin{2}; rad_probe = varargin{3};
        max_pos = varargin{4};
    case 5
        max_sz = varargin{1}; min_sz = varargin{2}; rad_probe = varargin{3};
        max_pos = varargin{4}; rot_angle = varargin{5};
    case 6
        max_sz = varargin{1}; min_sz = varargin{2}; rad_probe = varargin{3};
        max_pos = varargin{4}; rot_angle = varargin{5}; pos_angle = varargin{6};
    case 7
        max_sz = varargin{1}; min_sz = varargin{2}; rad_probe = varargin{3};
        max_pos = varargin{4}; rot_angle = varargin{5}; pos_angle = varargin{6};
        cond_r = varargin{7};
end
% Make cell arrays
opts_in = {max_sz, min_sz, rad_probe, max_pos, rot_angle, pos_angle, cond_r};
opts_out = cell( length(opts_in), n_objs);
for i = 1:1:length( opts_in );
    tmp = opts_in{i};
    if ~iscell(tmp)
        tmp = repmat( {tmp},1, n_objs);
    elseif length(tmp)==1
        tmp = repmat( tmp,1, n_objs); % is already a cell array
    elseif length(tmp)~=n_objs
        error('Parameters not understood');
    end
    for j=1:1:n_objs
        opts_out{i,j} = tmp{j};
    end
end
argout = opts_out;
end

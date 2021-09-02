function [ conds, info_obj ] = createMultiTarg( n_imgs, n_objs, overlap, mdl, bk, fcn, varargin )
%CREATEMULTITARG Creates multiple targets in one single FEM
%   Use this file to test a specific case in which there are more than one
%   object to detect
%	You can use a cell array to specify different functions to draw the
%	targets
%	Use cell arrays to get different properties for each target
%	@input
%		overlap: Whether or not the targets may overlap
%		fcn: The function to draw the targets
%			Inputs of this function will be: 1 (# images), fmdl, bk, varargin
%			default: @createEllipseGeneral
%		varargin: The rest of the parameters to generate the targets.
%		Please refer to the documentation of function fcn for information
%		on the input format. Use a cell array to specify different
%		properties for different targets
%	@output
%		conds: The list of conductivities of each element in each image
%		info_obj: array containing informations of random objects
%
%	Works for 2D applications. TODO: 3D
%
%	(C) 2015/06/15 Sebastien Martin


debug = false;

try	mdl2=mdl.fwd_model; catch; mdl2=mdl; end; mdl=mdl2;
if isfield(mdl,'coarse2fine'); mdl = rmfield(mdl,'coarse2fine'); end
if nargin<5; fcn = @createEllipseGeneral; end
[fcn, targ_props] = parse_varargin(n_objs, fcn, varargin{:}); % Target properties

n_elems = size(mdl.elems,1);
is2D = (size(mdl.nodes)==2);

conds = zeros(n_elems,n_imgs);

for k=1:1:n_imgs
    c_cond = zeros(n_elems,1);
    l=1; prev_objs = [];
    while (l<=n_objs)
        [new_cond, new_obj{l}] = fcn{l}(1,mdl,bk,targ_props{:,l});
        elems_obj = find( new_cond(:)~=bk ); % Get elements representing the object
        if ~overlap && sum( ismember(elems_obj,prev_objs)>0 ) % Check overlapping
            continue;
        else
            prev_objs = [prev_objs; elems_obj];
        end
        new_cond = new_cond-bk;
        % Here, it depends on how do you want to deal with overlapping
        % 		c_cond( elems_obj ) = c_cond(elems_obj)+new_cond(elems_obj);
        c_cond( elems_obj ) = new_cond(elems_obj);
        info_obj{k}(:,l) = new_obj{l};
        
        if debug
            figure; show_fem( mk_image(mdl, c_cond) );
        end
        l=l+1;
    end
    conds(:,k) = c_cond + bk;
    if debug
        figure; show_fem( mk_image(mdl, conds(:,k)), 1 ); title('Generated image');
    end
end

end

% Different properties for each object
% Errors are commented so we can easily change the number of targets we want
function [fcn, argout] = parse_varargin( n_objs, fcn, varargin )
% Get functions
if ~iscell(fcn); fcn = {fcn}; end
if length(fcn)==1; fcn = repmat(fcn,n_objs,1); end
% Make cell arrays
argout = cell( length(varargin), n_objs);
for i = 1:1:length( varargin );
    tmp = varargin{i};
    if ~iscell(tmp)
        tmp = repmat( {tmp},1, n_objs);
    elseif length(tmp)==1
        tmp = repmat( tmp,1, n_objs); % is already a cell array
    end
    for j=1:1:n_objs
        argout{i,j} = tmp{j};
    end
end
end

function [ n_mdl ] = shape_deformation( mdl, new_x, new_y )
%SHAPE_DEFORMATION Modify the shape of the mesh
%   In case of real applications, the FEM is not perfectly circular. This
%   file slightly updates the global shape to make it oval
%
%	WARNING: Generates oval shapes by extending the x or y axis
%	TO DO: 3D meshes, add z-axis
%
%	@input
%		mdl: a forward or inverse model
%		new_x: The new maximal X-axis coordinate
%		new_y: The new maximal Y-axis coordinate
%	@output
%		n_mdl: The deformed forward model
%
%	2014, Sebastien Martin

debug = false;
if isfield(mdl,'fwd_model'); mdl=mdl.fwd_model; end
n_mdl = mdl;

cur_x = max(mdl.nodes(:,1));
cur_y = max(mdl.nodes(:,2));
adj_x = new_x/cur_x;
adj_y = new_y/cur_y;
n_mdl.nodes(:,1) = n_mdl.nodes(:,1).*adj_x;
n_mdl.nodes(:,2) = n_mdl.nodes(:,2).*adj_y;

if debug
    figure;
    show_fem(n_mdl);
    title('Shape adjusted');
end
end
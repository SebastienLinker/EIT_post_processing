function [ fmdl ] = redo_elec_pos( fmdl, elec_pos, elec_shape, elec_obj, pb_rad)
%REDO_ELEC_POS Sometimes EIDORS gets it wrong
%   Idea: Check electrode position, electrode by electrode, and see if it
%   makes sense
%	Electrode is on 100+ nodes: Does not make sense
%	Electrode is not at the boundary of the FEM: Does not make sense (in our application)

debug = false;
n_elecs = length(fmdl.electrode);
if n_elecs~=length(elec_obj)
    error('Inconsistent sizes of electrodes');
end
if nargin < 5; pb_rad=1; end

% Find boundary nodes
bnd_nodes = find_boundary(fmdl.elems); bnd_nodes = unique(bnd_nodes(:));
eucl = sqrt( fmdl.nodes(:,1).^2 + fmdl.nodes(:,2).^2);
probe_nodes = find( (eucl<=pb_rad+10*eps) & (eucl>=pb_rad-100*eps));

for k=1:1:n_elecs
    curr_nodes = fmdl.electrode(k).nodes;
    if ~(all(ismember(curr_nodes,bnd_nodes)) && all(ismember(curr_nodes,probe_nodes)))
        % Here we need to adjust this electrode
        new_centre = elec_pos(k,[1:3]);
        new_nodes = find_elec_centre_internal(fmdl.boundary, fmdl.nodes,...
            fmdl.boundary_numbers, new_centre, probe_nodes, elec_shape);
        fmdl.electrode(k).nodes = new_nodes';
        if debug;
            figure; show_fem(fmdl); disp(['Electrode ',int2str(k),' has been repositioned']);
        end
        curr_nodes = new_nodes;
    end
    if length(curr_nodes)>50
        eidors_msg(['Electrode ',int2str(k),' is probably badly positioned'],1);
    end
end

end

% Inspired from ng_tank_find_elec, and simplified
% Don't understand the code, no comment, then just rewrite everything
function [nodes_elec] = find_elec_centre_internal(srf,vtx,fc,centre, pb_nodes, elec_shape)
p_nodes = vtx(pb_nodes,:); % Possible nodes to which the electrode is attached
if elec_shape(2)>0
    % Considers square electrodes
    dist_w = sqrt( (p_nodes(:,1) - centre(1,1)).^2 + (p_nodes(:,2) - centre(1,2)).^2 );
    nodes_w = (dist_w <= elec_shape(1)/2); % Width is (more than?) ok here
    dist_x = abs(p_nodes(:,1) - centre(1,1));
    nodes_x = (dist_x <= (10*eps +elec_shape(1)/2)); % MAYBE X-axis is ok here
    dist_y = abs(p_nodes(:,2) - centre(1,2));
    nodes_y = (dist_y <= (10*eps +elec_shape(1)/2)); % MAYBE Y-axis is ok here
    xy_nodes = pb_nodes(nodes_x & nodes_y);
    
    dist_z = abs(p_nodes(:,3) - centre(1,3));
    nodes_z = (dist_z <= (10*eps +elec_shape(2)/2)); % Height is ok here
    nodes_elec = pb_nodes(nodes_x & nodes_y & nodes_z);
else
    error(['You probably have selected circular electrodes and one',...
        'of them seems to be poorly placed']);
end

end
function [ area_elems ] = calc_elements_area( mdl )
%CALC_ELEMENTS_AREA Calculate the area or volume of each element in the model
%   Input:
%		mdl: the forward model
%	Output:
%		area_elems: row vector containing the area (or volume) of each element
%
%	Am I doing the same thing as 'get_elem_volume' ?

area_elems = eidors_obj('get-cache', mdl, 'area_elem');
if ~isempty(area_elems); return; end

n_elems = size(mdl.elems,1);

if (size(mdl.elems,2)==3) %2D model
    nodes_vec = zeros( n_elems, size(mdl.elems,2)*2);
    nodes_vec(:,[1,2]) = mdl.nodes(mdl.elems(:,1),:);
    nodes_vec(:,[3,4]) = mdl.nodes(mdl.elems(:,2),:);
    nodes_vec(:,[5,6]) = mdl.nodes(mdl.elems(:,3),:);
    area_elems = triangleArea(nodes_vec(:,[1 2]), nodes_vec(:,[3 4]), nodes_vec(:,[5 6]));
elseif (size(mdl.elems,2)==4) %3D model
    nodes_vec = zeros( n_elems, size(mdl.elems,2)*3);
    nodes_vec(:,[1,2,3]) = mdl.nodes(mdl.elems(:,1),:);
    nodes_vec(:,[4,5,6]) = mdl.nodes(mdl.elems(:,2),:);
    nodes_vec(:,[7,8,9]) = mdl.nodes(mdl.elems(:,3),:);
    nodes_vec(:,[10,11,12]) = mdl.nodes(mdl.elems(:,4),:);
    for k=1:1:n_elems
        area_elems(k) = tetrahedronVolume([nodes_vec(k,[1 2 3]); nodes_vec(k,[4 5 6]); nodes_vec(k,[7 8 9]); nodes_vec(k,[10 11 12])] );
    end
end

eidors_obj('set-cache', mdl, 'area_elem', area_elems);

end
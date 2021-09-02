function [ isopen ] = isOpenDomain( mdl )
%ISOPENDOMAIN Determines whether the model is open (or unbounded) domain or
%close domain EIT
%   Basically checks the position of each electrode and checks whether
%   there are nodes farther than the electrodes
%	May not work very well if the probe is on the side of the domain
%	(though this would be a strange configuration)
%
%	(C) 2015/08/19 Sebastien Martin

if isfield(mdl,'fwd_model'); mdl = mdl.fwd_model; end


if isfield(mdl,'parameters') && isfield(mdl.parameters,'isOpenDomain')
    isopen = mdl.parameters.isOpenDomain;
    return;
end

nodes_elecs = [mdl.electrode(:).nodes];
pos_elec = mdl.nodes(nodes_elecs,:);
dist_elecs = sqrt( pos_elec(:,1).^2 + pos_elec(:,2).^2 );

max_dist_node = max( sqrt( mdl.nodes(:,1).^2 + mdl.nodes(:,2).^2 ) );
isopen = max_dist_node*0.9 > max(dist_elecs);

end
function [ conds, info_obj ] = createConvex2Dmulti_obj( n_objs, n_imgs, mdl, bk, varargin )
%CREATECONVEX2DMULTI_OBJ Creates multiple convex shapes in one single FEM
%   Use this file to test a specific case in which there are more than one
%   object to detect
%	@input
%	@output conds: The list of conductivities of each element in each image
%		info_obj: array containing size, position(X,Y) and conductivity of
%		random objects (having the highest and lowest conductivities): Total size= n_imgs*8

debug = false;

try	mdl2=mdl.fwd_model; catch; mdl2=mdl; end; mdl=mdl2;
n_elems = size(mdl.elems,1);

conds = zeros(n_elems,n_imgs);

for k=1:1:n_imgs
    c_cond = zeros(n_elems,1);
    l=1; prev_objs = [];
    while (l<=n_objs)
        [new_cond, new_obj{l}] = createConvex(1,mdl,bk,varargin{:});
        elems_obj = find( new_cond(:)~=bk ); % Get elements representing the object
        prev_objs = [prev_objs; elems_obj];
        new_cond = new_cond-bk;
        % Here, it depends on how do you want to deal with overlapping
        % 		c_cond( elems_obj ) = c_cond(elems_obj)+new_cond(elems_obj);
        c_cond( elems_obj ) = new_cond(elems_obj);
        
        if(l>1)
            if (new_obj{l}(4)>info_obj(4,k))
                info_obj([1:4],k) = new_obj{l};
            end
            if (new_obj{l}(4)<info_obj(8,k))
                info_obj([5:8],k) = new_obj{l};
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


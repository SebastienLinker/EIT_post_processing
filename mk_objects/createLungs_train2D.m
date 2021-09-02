function [ conds, info_obj, all_info ] = createLungs_train2D( n_imgs, mdl, bk, cond_r, n_objs )
%CREATELUNGS_TRAIN2D Creates multiple convex shapes, imitating EIT of lungs
%	This is made for training purpose, so the size is variable
%   Code based on function createConvex2Dmulti_obj
%	@input
%		n_imgs
%		mdl
%		bk: background conductivity
%		cond_r1: Conductivity range for the lungs, usually [inspiration expiration]
%	@output conds: The list of conductivities of each element in each image
%		info_obj: array containing conductivity of the lungs
%		all_info: Array with size (x,y) CoG (x,y) and cond. of each lung
%

debug = false;
if debug; warning('Debug mode activated'); end

if ~exist('cond_r2','var'); cond_r2=cond_r; end;
if ~exist('n_objs','var'); n_objs=2; end; % 2 lungs + one blood clot

try	mdl2=mdl.fwd_model; catch; mdl2=mdl; end; mdl=mdl2;
if isfield(mdl,'coarse2fine'); mdl = rmfield(mdl,'coarse2fine'); end
n_elems = size(mdl.elems,1);

conds = zeros(n_elems,n_imgs);
lst_cond = rand(n_imgs,1);

% TODO: caller should specify this separately in passed arguments
cond_blood = rand(n_imgs,1) * (cond_r(1,2)-cond_r(1,1))+cond_r(1,1);
if size(cond_r,1)>=2
    cond_r = cond_r(2:1:end,:);
end

tot_range = sum(cond_r(:,2)-cond_r(:,1));
cond_th = [0 ((cond_r(:,2)-cond_r(:,1))/tot_range)'];
for k=2:1:size(cond_r,1)+1
    cond_th(k) = cond_th(k)+cond_th(k-1);
    tmp_elems = find( (lst_cond>=cond_th(k-1)) & (lst_cond<=cond_th(k)) );
    tmp_inhomo(tmp_elems) = lst_cond(tmp_elems)* ...
        (cond_r(k-1,2)-cond_r(k-1,1))+cond_r(k-1,1) - bk; % xxx*(max-min)+min
end
lst_cond = tmp_inhomo' + bk;

for k=1:1:n_imgs
    c_cond = zeros(n_elems,1);
    l=1; prev_objs = [];
    while (l<=n_objs-1)
        if (l==1) % Lungs
            [new_cond, new_obj{l}, CoGL{k}] = createLungs(mdl,bk, lst_cond(k));
        else % Circular objects (blood)
            [new_cond, tmp] = createConvex(1,mdl,bk, 0.2,0.1, 0,0.5, 0, ...
                360, [cond_blood(k), cond_blood(k)]);
            new_obj{l} = [tmp(1) tmp];
        end
        
        elems_obj = find( new_cond(:)~=bk ); % Get elements representing the object
        if ((l>1) && sum( ~ismember(elems_obj,prev_objs)>0 )) % Objects should overlap each other
            continue;
        else
            prev_objs = [prev_objs; elems_obj];
        end
        new_cond = new_cond-bk;
        c_cond( elems_obj ) = new_cond(elems_obj);
        
        if(l>1)
            if (new_obj{l}(5)>info_obj(5,k)) % Higher conductivity
                info_obj(1:5,k) = new_obj{l};
            end
            if (new_obj{l}(5)<info_obj(10,k)) % Higher conductivity
                info_obj(6:10,k) = new_obj{l};
            end
        else
            info_obj(:,k) = [new_obj{1}, new_obj{1}];
        end
        
        if debug
            figure; show_fem( mk_image(mdl, c_cond),1 );
        end
        l=l+1;
    end
    
    conds(:,k) = c_cond + bk;
    if debug
        figure; show_fem( mk_image(mdl, conds(:,k)), 1 );
    end
end

all_info = CoGL;
info_obj = info_obj([5 10],:)';
end

%% Create two lungs
% Algorithm based on createConvex and createEllipse
function [ img_data, info_obj, CoGs ] = createLungs(mdl, bk, cond)
n_elems = size(mdl.elems,1);

params = rand(1,3);
X_sz = max(mdl.nodes(:,1)); Y_sz = max(mdl.nodes(:,2));
szX = (params(1)*0.2+0.25)*X_sz; szY = (params(1)*0.3+0.35)*Y_sz;
posX(1) = (params(2)*0.2-0.1 -0.5)*X_sz; % Left lung
posX(2) = (params(2)*0.2-0.1 +0.5)*Y_sz; % Right lung
posY = 0;

elems = bk * ones(n_elems,1);
info_obj = [szX, szY, posX(1), posY, cond]; % Conductivity should be the 5th one
CoGs = [posX(1), posY, szX; posX(2), posY, szY]';

% Thresholding algorithm
inhomo = mk_image( mdl, elems );
k=1; img_data = inhomo.elem_data;
while (k<=2)
    cur_elems{k} = zeros(1, n_elems);
    select_fcn = inline(['((x-',num2str(posX(k)),')/',num2str(szX),').^2',...
        '+((y-',num2str(posY),')/',num2str(szY),').^2','<1'],'x','y','z');
    cur_elems{k}(1,:) = bk - cond * elem_select( inhomo.fwd_model, select_fcn );
    threshold = bk - cond*0.1;
    % variable conductivity range
    elems_obj = cur_elems{k}(1,:)<threshold;
    elems_no_obj = cur_elems{k}(1,:)>threshold;
    cur_elems{k}(1, elems_obj) = cond;
    cur_elems{k}(1, elems_no_obj) = bk;
    if( min(cur_elems{k})==max(cur_elems{k}) ) % homogeneous image? (due to conductivity)
        warning('It seems that the algorithm failed to generate the lung');
    end
    
    % Check overlapping (should not be, but just in case
    if ((k==2) && sum( ismember( find(cur_elems{1}~=bk), find(cur_elems{2}~=bk) )>0))
        [img_data, info_obj, CoGs] = createLungs(mdl, bk, cond);
        return;
    end
    img_data(elems_obj) = cur_elems{k}(1, elems_obj==1 )';
    k=k+1;
end

end






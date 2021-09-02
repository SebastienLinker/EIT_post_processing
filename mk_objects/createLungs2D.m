function [ conds, info_obj ] = createLungs2D( n_imgs, mdl, bk, cond_r1, cond_r2 )
%CREATELUNGS2D Creates multiple convex shapes, imitating EIT of lungs
%   Based on function createConvex2Dmulti_obj
%	@input
%		n_imgs
%		mdl
%		bk: background conductivity
%		cond_r1: Conductivity range for the lungs, usually [inspiration expiration]
%	@output conds: The list of conductivities of each element in each image
%		info_obj: array containing conductivity of the lungs
%

debug = false;
if debug; warning('Debug mode activated'); end

if ~exist('cond_r2','var')
    cond_r2=cond_r1;
end

try	mdl2=mdl.fwd_model; catch; mdl2=mdl; end; mdl=mdl2;
n_elems = size(mdl.elems,1);
n_objs=2; % 2 lungs + one blood clot

conds = zeros(n_elems,n_imgs);
lst_cond = rand(n_imgs,1);
tot_range = sum(cond_r1(:,2)-cond_r1(:,1));
cond_th = [0 ((cond_r1(:,2)-cond_r1(:,1))/tot_range)'];
for k=2:1:size(cond_r1,1)+1
    cond_th(k) = cond_th(k)+cond_th(k-1);
    tmp_elems = find( (lst_cond>=cond_th(k-1)) & (lst_cond<=cond_th(k)) );
    tmp_inhomo(tmp_elems) = lst_cond(tmp_elems)* ...
        (cond_r1(k-1,2)-cond_r1(k-1,1))+cond_r1(k-1,1) - bk; % xxx*(max-min)+min
end
lst_cond = tmp_inhomo' + bk;

for k=1:1:n_imgs
    c_cond = zeros(n_elems,1);
    l=1; prev_objs = [];
    while (l<=n_objs)
        if (l<=2) % Lungs
            if (l==1) %Left lungs
                [new_cond, new_obj{l}] = createLung(mdl,bk, 'left', lst_cond(k));
            else %Right lungs
                [new_cond, new_obj{l}] = createLung(mdl,bk, 'right', lst_cond(k));
            end
        else % Circular objects (blood)
            [new_cond, tmp] = createConvex(1,mdl,bk, 0.2,0.05, 0,0.5, 0,360, cond_r2);
            new_obj{l} = [tmp(1) tmp];
        end
        
        elems_obj = find( new_cond(:)~=bk ); % Get elements representing the object
        if ((l>2) && sum( ~ismember(elems_obj,prev_objs)>0 )) % Objects should overlap each other
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
            info_obj(:,k) = [new_obj{l}, new_obj{l}];
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

info_obj = info_obj([5 10],:)';
end

%% Create one lung
% Algorithm based on createConvex and createEllipse
function [ img_data, info_obj ] = createLung(mdl, bk, side, cond)
n_elems = size(mdl.elems,1);
[szX, szY, posX, posY] = regenLung( side );
elems = bk * ones(n_elems,1);
info_obj = [szX, szY, posX, posY, cond];

% Thresholding algorithm
inhomo = mk_image( mdl, elems );
elems1 = zeros(1, n_elems);
% 	for k = 1:1:n_images
select_fcn = inline(['((x-',num2str(posX),')/',num2str(szX),').^2',...
    '+((y-',num2str(posY),')/',num2str(szY),').^2',...
    '<1'],'x','y','z');

elems1(1,:) = bk - cond * elem_select( inhomo.fwd_model, select_fcn );
threshold = bk - cond*0.1;
% variable conductivity range
elems_obj = find(elems1(1,:)<threshold);
elems_no_obj = find(elems1(1,:)>threshold);
elems1(1, elems_obj) = cond; % bkgnd - inhomo_conduct(k);
elems1(1, elems_no_obj) = bk;

img_data = elems1(1,:)';
% 	end
% homogeneous image? (due to conductivity)
if( min(img_data)==max(img_data) )
    warning('It seems that the algorithm failed to generate the lung');
end
end

function [szX, szY, posX, posY] = regenLung( side )
szX=0.3; szY=0.6; posY=0;
if strcmp(side,'left'); posX =-0.5; else posX = 0.5;	end
end






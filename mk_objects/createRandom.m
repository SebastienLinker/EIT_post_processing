function [ img_data, info_obj, targs ] = createRandom( n_imgs, mdl, varargin )
%CREATERANDOM Generate targets, random shape
%   This function generates n_imgs EIDORS images, each of them containing a
%   specific target
%	Targets can be any of the target generated by other functions called
%	'createSomething'
%
%	Inputs
%		n_imgs: Number of images you need
%		mdl: Inverse or forward model
%		varargin: Additional arguments, sent to target creation functions
%		as is
%	Outputs
%		img_data: list of conductivities
%		info_obj: A N*n_imgs array containing the information on each
%		conductivity distribution
%		eqn: A cell array containing the equations of each target (inline
%		function)
%		targs: The list of each target (cell array)
%
%	(C) 2015/08/10 Sebastien Martin

debug = false;

try mdl = mdl.fwd_model; catch; end
is2D = (size(mdl.nodes,2)==2);
n_elems = size(mdl.elems,1);

[lst_fcn, N_fcn] = getListOfExistingFunctions(is2D);
assigned_num = ceil(rand(n_imgs,1).*N_fcn);
targs = lst_fcn(assigned_num);

% Allocate memory
img_data = zeros(n_elems, n_imgs);
info_obj = zeros(n_imgs, 9);
eqn = cell(n_imgs,1);
n_infos_max = 0;

for k=1:1:N_fcn
    c_imgs = (assigned_num == k); % Current images
    tmp_n_imgs = sum(c_imgs);
    if tmp_n_imgs==0; continue; end; %No image with this shape
    
    [tmp_conds, tmp_info] = lst_fcn{k}( tmp_n_imgs, mdl, varargin{:}); % Generate targets
    
    if debug
        disp(['Generated ',int2str(tmp_n_imgs),' images with ',func2str(lst_fcn{k})]);
        disp_debug(mdl,tmp_conds);
    end
    
    img_data(:, c_imgs) = tmp_conds;
    n_infos = size(tmp_info,2);
    n_infos_max = max(n_infos_max, n_infos);
    if n_infos<9
        tmp_info = [tmp_info, zeros(tmp_n_imgs, 9-n_infos)];
    end
    info_obj( c_imgs,:) = tmp_info;
end

% Reduce size of info_obj array
zero_column = all(info_obj==0,1);
info_obj = info_obj(:,~zero_column);

end

function [ list, N ] = getListOfExistingFunctions(is2D)
if is2D
    list = {
        @createEllipseGeneral;
        @createTriangle;
        @createSquare;
        };
else
    list = {
        @createCylinder;
        @createEllipseGeneral;
        @createTriangle;
        @createSquare;
        };
end
N = length(list);
end

function [ list ] = getListOfRealisticShapes(is2D)
if is2D
    list = {@createLungs_train2D};
else
    list = [];
end
end

function [] = disp_debug(mdl, conds)
n_imgs = size(conds,2);
for k=1:1:n_imgs
    img = mk_image(mdl, conds(:,k));
    figure; show_fem(img,1); title('Debug mode');
end
end
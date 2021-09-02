function [ GREIT, AR, PE, RES, SD, RNG ] = cmp_multi_3D_GREIT( imgs, orig, f_high, varargin )
%CMP_MULTI_3D_GREIT Compare GREIT errors of a set of images instead of a
%single image
%	Inputs
%		imgs: set of images, can be either matrix or cell array
%		orig: The original image
%		f_high: Force the highest peak to be positive, may reverse
%			everything in unbounded domain EIT (default true)
%		pts: number (or vector) of points
%	Outputs
%		GREIT: The outputs from eval_GREIT_fig_merit or cmp_GREIT_fig_merit
%
%	2015/10/26 Sebastien Martin

debug = false;
if nargin<=2
    f_high = true;
end
n_imgs = numel(imgs);
sz_imgs = size(imgs);

% Verify input data
if ismatrix(imgs)
    imgs = num2cell(imgs);
end

if size(orig) == [1 1]
    orig = repmat(orig,1,n_imgs);
end
orig = mat2cell(orig,1,ones(n_imgs,1));

% Adjust the conductivity range
imgs = cellfun(@(x,y,z) preprocess_greit(x, f_high, varargin{:}), {imgs{:}}, 'UniformOutput',false);
orig = cellfun(@(x,y,z) preprocess_greit(x, f_high, varargin{:}), {orig{:}}, 'UniformOutput',false);
% Actually get the GREIT errors
GREIT = cell2mat( cellfun(@(x,y) cmp_GREIT_3D(x,y), imgs, orig, 'UniformOutput',false));

% Extract, AR, PE, RES, SD, and RNG separately
ext_errs = @(x) reshape(GREIT(x,:),sz_imgs);
AR = ext_errs(1); PE = ext_errs(2); RES =ext_errs(3); SD = ext_errs(4); RNG = ext_errs(5);

if debug
    disp_dbg(imgs);
    disp_dbg(orig);
end

end

function img = preprocess_greit(img, f_high, pts, area)
% Make sure the background is 0 and the target ocnductivity is above the background
img_data = get_img_data(img);
if ~ischar(f_high) && (f_high==true)
    img_data = img_data - median(img_data);
    img.elem_data = img_data;
    img = force_direction(img,'higher');
elseif ~ischar(f_high)
    % 		img.elem_data = img.elem_data-median(img.elem_data);
    area_vec = abs(calc_elements_area(img.fwd_model));
    if ~all(size(area_vec)==size(img_data))
        area_vec = area_vec';
        if ~all(size(area_vec)==size(img_data))
            error('Please verify your image data'); end
    end
    mean_vec = img_data .* area_vec;
    img_data = img_data - (mean(mean_vec)/mean(area_vec));
    img.elem_data = img_data;
elseif strcmp(f_high,'no_norm')
    % No normalization at all, leave everythin as is
else
    error('Unknown parameter, f_high');
end

% Specify number of points
if nargin>3
    switch length(pts)
        case 1
            img.fwd_model.mdl_slice_mapper.npx = pts; img.fwd_model.mdl_slice_mapper.npy = pts;
        case 2
            img.fwd_model.mdl_slice_mapper.npx = pts(1); img.fwd_model.mdl_slice_mapper.npy = pts(2);
        otherwise
            if min(size(pts))==1
                img.fwd_model.mdl_slice_mapper.x_pts = pts; img.fwd_model.mdl_slice_mapper.y_pts = pts;
            else
                img.fwd_model.mdl_slice_mapper.x_pts=pts(1,:);img.fwd_model.mdl_slice_mapper.y_pts=pts(2,:);
            end
    end
end
% If a specific area is defined
if nargin > 4
    X = linspace(area(1), area(2), img.fwd_model.mdl_slice_mapper.npx);
    Y = linspace(area(3), area(4), img.fwd_model.mdl_slice_mapper.npy);
    img.fwd_model.mdl_slice_mapper.x_pts = X;
    img.fwd_model.mdl_slice_mapper.y_pts = Y;
    img.fwd_model.mdl_slice_mapper = rmfield(img.fwd_model.mdl_slice_mapper,{'npx','npy'});
end
end

function [] = disp_dbg(imgs)
for k = 1:1:length(imgs)
    img = imgs{k};
    figure; show_fem(img,1); axis off;
end
end
function [ GREIT, AR, PE, RES, SD, RNG ] = cmp_multi_GREIT( imgs, orig, xyzr, f_high, varargin )
%CMP_MULTI_GREIT Compare GREIT errors of a set of images instead of a
%single image
%   While it is possible to do so directly within cmp_GREIT_fig_merit,
%   some cases, especially 3D, may lead to a crash
%	It is adviced to use this wrapper instead
%	Inputs
%		imgs: set of images, can be either matrix or cell array
%		orig: The original image
%		Z: Coordinates of the CoGs of the target: a 4*1 or 4*N matrix
%			For 2D applications, a 3*1 or 3*N matrix is ok
%			If it is an image, then cmp_GREIT_fig_merit is called and this
%			image is considered as reference image
%		f_high: Force the highest peak to be positive, may reverse
%			everything in unbounded domain EIT (default true)
%		pts: number (or vector) of points
%	Outputs
%		GREIT: The outputs from eval_GREIT_fig_merit or cmp_GREIT_fig_merit
%
%	2015/07/21 Sebastien Martin

debug = false;
if nargin<=3; f_high = true; end
n_imgs = numel(imgs); sz_imgs = size(imgs);

% Verify input data
if ismatrix(imgs); imgs = num2cell(imgs); end
if isstruct(orig); use_depr_call = true; else use_depr_call = false; end

if use_depr_call
    if size(orig) == [1 1]; orig = repmat(orig,1,n_imgs); end
    orig = mat2cell(orig,1,ones(n_imgs,1));
    if size(xyzr,2)==1; xyzr = repmat(xyzr,1,n_imgs); end
    if size(xyzr,1)==3; xyzr = [xyzr([1 2],:) zeros(1,n_imgs) xyzr(3,:)]; end
    xyzr = mat2cell(xyzr,4,ones(n_imgs,1));
else
    t_type = xyzr{2}; xyzr = xyzr{1};
end

% Adjust the conductivity range
if use_depr_call
    imgs = cellfun(@(x,y,z) preprocess_greit(x,y(3), f_high, varargin{:}), {imgs{:}}, xyzr, 'UniformOutput',false);
    orig = cellfun(@(x,y,z) preprocess_greit(x,y(3), f_high, varargin{:}), {orig{:}}, xyzr, 'UniformOutput',false);
    % Actually get the GREIT errors
    GREIT = cell2mat( cellfun(@(x,y) cmp_GREIT_fig_merit(x,y), imgs, orig, 'UniformOutput',false));
else
    imgs = cellfun(@(x) preprocess_greit(x,xyzr(3), f_high, varargin{:}), {imgs{:}}, 'UniformOutput',false);
    % Actually get the GREIT errors
    GREIT = cell2mat( cellfun(@(x) cmp_GREIT_fig_merit(x, xyzr, t_type), imgs, 'UniformOutput',false));
end

% Extract, AR, PE, RES, SD, and RNG separately
ext_errs = @(x) reshape(GREIT(x,:),sz_imgs);
AR = ext_errs(1); PE = ext_errs(2); RES =ext_errs(3); SD = ext_errs(4); RNG = ext_errs(5);

if debug
    disp_dbg(imgs, xyzr);
    if use_depr_call; disp_dbg(orig, xyzr); end
end

end

function img = preprocess_greit(img, Z, f_high, pts, area)
% Make sure the background is 0 and the target conductivity is above the background
if ~ischar(f_high) && (f_high==true)
    img.elem_data = img.elem_data - median(img.elem_data);
    img = force_direction(img,'higher');
elseif ~ischar(f_high)
    % 		img.elem_data = img.elem_data-median(img.elem_data);
    area_vec = abs(calc_elements_area(img.fwd_model));
    if ~all(size(area_vec)==size(img.elem_data))
        area_vec = area_vec';
        if ~all(size(area_vec)==size(img.elem_data))
            error('Please verify your image data'); end
    end
    mean_vec = img.elem_data .* area_vec;
    img.elem_data = img.elem_data- (mean(mean_vec)/mean(area_vec));
elseif strcmp(f_high,'no_norm')
    % No normalization at all, leave everythin as is
else
    error('Unknown parameter, f_high');
end
% If 3D model, we need to specify the Z-axis
if size(img.fwd_model.nodes,2)==3
    img.calc_slices.levels = [Inf Inf Z];
end

% Specify number of points
if nargin>3 && ~iscell(pts) && ~ischar(pts)
    switch length(pts)
        case 1
            img.fwd_model.mdl_slice_mapper.npx = pts; img.fwd_model.mdl_slice_mapper.npy = pts;
        case 2
            img.fwd_model.mdl_slice_mapper.npx = pts(1); img.fwd_model.mdl_slice_mapper.npy = pts(2);
        otherwise
            if min(size(pts))==1;
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

function [] = disp_dbg(imgs, xyzr)
for k = 1:1:length(imgs)
    img = imgs{k}; pos = xyzr{k}([1 2]); rad = xyzr{k}(4);
    slice = calc_slices(img);
    sz = size(slice);
    figure;
    subtightplot(1,2,1); show_slices(slice);
    subtightplot(1,2,2); show_fem(img,1); axis off;
    suptitle([num2str(sz(1)),'*',num2str(sz(2)),' image, Target size: ',num2str(rad),...
        ' Location: ',num2str(pos(1))]);
end
end
function [ img_err, n_err, avg_err_dist, err_area_pct, th ] = area_error_unbounded( targ, rec, pos, varargin )
%AREA_ERROR_UNBOUNDED Compute the area error for unbounded EIT images
%   In unbounded EIT, we define the area error as follow:
%	The image is reduced to keep only the elements between the target (CoG)
%	and the probe (CoG). A second ellipse is generated, similar to the
%	first one, but with larger radii. Error elements also belong to this
%	second ellipse
%	The conductivity of each element is compared to a threshold
%	The area error is the number of elements for which the conductivity in
%	the reconstruction is below the threshold, and above in the original,
%	or the opposite
%	Returns the images, the number of elements considered as error,
%				the average error in distance, and the percentage of error.
%	Input:
%		targ: target image (simulation, normalized), can be a vector
%		rec: reconstruction (after normalization)
%			If rec is a matrix, then one line of this matrix corresponds to
%			the same target image. Matrices can be used, for example with
%		pos: properties of the ellipse, the output of createEllipseGeneral (N*6) or (N*9) array
%		varargin: A pair 'name','value' of arguments. Name can take the
%			following values:
%			slice: Next argument is a 3*1 vector indicating the slice
%			chk_dir: Next argument is a boolean, decides to check the
%				direction or not. Default false
%	Output:
%		img_err: EIDORS img containing NaN where an error is detected
%		n_err: # of elements containing an error
%		avg_err_dist: area of error elements / perimeter of arc of the semiellipse
%		err_area_pct: pct of the semiellipse considered as error (area)
%		th: the threshold determined from automatic estimate
%
%	(C) Sebastien Martin 2014/11/24

debug = false;

[slice_, chk_dir] = process_args(varargin{:});

if ~isvector(targ) || (length(targ)~=size(pos,1))
    if ~isvector(targ)
        error('Input ''targ'' is not a vector');
    else
        error('Number of target images doesn''t match the number of ellipses');
    end
end
if size(targ,1)==1; targ = targ'; end

% Check matrix rec
if (size(rec,1)~=size(targ,1))
    if (size(rec,2)==size(targ,1))
        rec = rec';
    else error('Inconsistent matrix sizes');
    end
end

for k=1:1:length(targ)
    for m=1:1:size(rec,2)
        [img_err(k,m), n_err(k,m), avg_err_dist(k,m), err_area_pct(k,m), th(k,m)] = ...
            area_error_internals( targ(k), rec(k,m), pos(k,:), slice_, chk_dir );
    end
end

end

function [ img_err, n_err, avg_err_dist, err_area_pct, th ] = ...
    area_error_internals( targ, rec,pos, slice_, chk_dir )

switched = [];
if chk_dir.do
    [targ, switched] = force_direction(targ, 'lower', 0);
    rec(switched) = force_direction(rec(switched), 'lower', 0);
end

debug = evalin('caller','debug');
do_slice = slice_.do;
if do_slice
    slice = slice_.slice;
end
fmdl = targ(1).fwd_model;
is2D = size(fmdl.nodes,2)==2;
area_elems = abs(calc_elements_area(rec.fwd_model));
th = (max(targ.elem_data)-min(targ.elem_data))*0.5;
if do_slice
    th = 0.25;
end
th(switched) = -th(switched);

x_pts = 128;
y_pts = 128;
rec.calc_colours.npoints = x_pts;
targ.calc_colours.npoints = x_pts;
det_th = rec.elem_data;
if debug
    figure; imhist( rec.elem_data ); title('Reconstruction');
    figure; imhist( det_th(~isnan(det_th)) ); title('Converted to image');
end
if is2D
    th_ell = [2 2]; % 1*2 array used to create the new ellipse: new radii = old radii .* th_ell
else
    th_ell = [2 2 2];
end

if do_slice
    if isfield(rec,'calc_colours') && isfield(rec.calc_colours,'npoints')
        X_npts = rec.calc_colours.npoints;
        Y_npts = rec.calc_colours.npoints;
    else
        X_npts = calc_colours('npoints');
        Y_npts = calc_colours('npoints');
    end
    if isfield(rec.fwd_model,'mdl_slice_mapper') && isfield(rec.fwd_model.mdl_slice_mapper,'x_pts')
        X_minus = min(rec.fwd_model.mdl_slice_mapper.x_pts);
        X_plus = max(rec.fwd_model.mdl_slice_mapper.x_pts);
        X_npts = length(rec.fwd_model.mdl_slice_mapper.x_pts);
    else
        X_minus = min(fmdl.nodes(:,1)); X_plus = max(fmdl.nodes(:,1));
    end
    if isfield(rec.fwd_model,'mdl_slice_mapper') && isfield(rec.fwd_model.mdl_slice_mapper,'y_pts')
        Y_minus = min(rec.fwd_model.mdl_slice_mapper.y_pts);
        Y_plus = max(rec.fwd_model.mdl_slice_mapper.y_pts);
        Y_npts = length(rec.fwd_model.mdl_slice_mapper.y_pts);
    else
        Y_minus = min(fmdl.nodes(:,2)); Y_plus = max(fmdl.nodes(:,2));
    end
    areaImage = (X_plus-X_minus)*(Y_plus-Y_minus);
    area1px = areaImage / (X_npts*Y_npts);
end

% Reduce image
sel_fcn = {};
if is2D % 2D mesh
    ell_c = pos(:,[4 5]); % centers
    ell_sz = pos(:,[1 2]); % To do: Consider the possible rotation
else
    ell_c = pos(:,[6 7 8]);
    ell_sz = pos(:,[1 2 3]);
end
eidors_msg('@@@ Warning: assumes the ellipse to be centered and aligned on x-axis',3);
% sel_fcn{end+1} = inline('x > 0','x','y','z');
% sel_fcn{end+1} = inline(['x < ',num2str(ell_c(1))],'x','y','z');
sel_fcn{end+1} = inline(['x > ',num2str(ell_c(1)-ell_sz(1)*th_ell(1))],'x','y','z');
sel_fcn{end+1} = inline(['x < ',num2str(ell_c(1)+ell_sz(1)*th_ell(1))],'x','y','z');
sel_fcn{end+1} = inline(['y > ',num2str(ell_c(2)-ell_sz(2)*th_ell(2))],'x','y','z');
sel_fcn{end+1} = inline(['y < ',num2str(ell_c(2)+ell_sz(2)*th_ell(2))],'x','y','z');
if ~is2D
    % Z-axis
end
memb_frac = elem_select( fmdl, sel_fcn);
new_sz = ell_sz.*th_ell;
if is2D
    new_ell = ~(createEllipseGeneral(1,fmdl,2,new_sz,new_sz,ell_c(1),ell_c(1),[pos(3),pos(3)],[0,0])==2);
else
    new_ell = ~(createEllipseGeneral(1,fmdl,2,new_sz,new_sz,ell_c([1 3]),ell_c([1 3]),...
        [pos(4),pos(4)],[0,0,0])==2);
end
memb_frac = (memb_frac==1) & (new_ell==1);

if debug
    figure;
    show_fem(mk_image(fmdl,memb_frac),1);
    title('Error region');
end

% Reduced targets and reconstructions: we only keep elements within the
% error region
if (is2D || ~do_slice)
    red_targ = targ.elem_data(memb_frac);
    red_rec = rec.elem_data(memb_frac);
    if debug
        figure;
        subtightplot(1,2,1);
        show_fem(targ,1);
        title('Target');
        subtightplot(1,2,2);
        show_fem(rec,1);
        title('Reconstruction');
    end
else
    % For 3D open domain, get one slice of the picture and estimate the
    % error on this slice
    mdl_frac = mk_image(fmdl,memb_frac);
    mdl_frac.calc_slices.levels = slice;
    mdl_frac.calc_colours.npoints = x_pts;
    targ.calc_slices.levels = slice;
    rec.calc_slices.levels = slice;
    red_targ = calc_slices(targ);
    red_rec = calc_slices(rec);
    memb_frac = calc_slices(mdl_frac);
    red_targ(isnan(red_targ)) = 0; red_rec(isnan(red_rec)) = 0; memb_frac(isnan(memb_frac)) = 0;
    if debug
        figure;
        subtightplot(1,3,1); show_slices(targ); title('Target');
        subtightplot(1,3,2); show_slices(rec); title('Reconstruction');
        subtightplot(1,3,3); show_slices(mdl_frac); title('ROI');
    end
end

% Calculate error on reduced images
% Get objects in original and output images
orig_inhomo = (red_targ < th);
result_inhomo = (red_rec < th);

if debug && ~do_slice
    img_inho = mk_image(fmdl, 0); img_rec = mk_image(fmdl, 0);
    img_inho.elem_data(memb_frac) = double(orig_inhomo);
    img_rec.elem_data(memb_frac) = double(result_inhomo);
    img_inho.calc_colours = struct('ref_level',0,'clim',1);
    img_rec.calc_colours = struct('ref_level',0,'clim',1);
    figure;
    subtightplot(1,2,1); show_fem(img_inho,1); title('Inhomo region, target image');
    subtightplot(1,2,2); show_fem(img_rec,1); title('Inhomo region, reconstruction image');
end
% Get elements in A,B, and C, where B is the region where both original and output images show an
% inhomogeneity, A is the region where only the original image has an
% inhomogeneity, C is the region where only the output image has one, and
% in D, there is no inhomogeneity, A and C contain the errors
elems_B = (orig_inhomo & result_inhomo);
elems_A = (orig_inhomo & ~result_inhomo);
elems_C = (~orig_inhomo & result_inhomo);
elems_D = (~orig_inhomo & ~result_inhomo);

n_err = sum(sum(elems_A)+sum(elems_C)); % # of elements containing an error

% Get error elements
if (is2D || ~do_slice); lst_elems = find(memb_frac==1);
else lst_elems = (memb_frac==1);
end
err_elems = lst_elems([elems_A | elems_C]);

% Image showing error elements
if (is2D || ~do_slice)
    img_err = rec; img_err.elem_data(err_elems) = NaN;
    if debug; figure; show_fem(img_err,1); end
    area_err = sum(area_elems(err_elems));
    area_frac = sum(area_elems(memb_frac));
else
    img_err = 0; % Don't output in this case
    if debug; figure; imshow(err_elems); end
    area_err = sum(err_elems)*area1px;
    area_frac = area_err/numel(err_elems);
end

% Perimeter of the original ellipse (initial image)
if is2D
    ell_per = ellipsePerimeter(ell_sz);
elseif do_slice
    % 	ell_per = sum(orig_inhomo(:))/2; %Seems working like that
    perim = bwperim(orig_inhomo);
    sz_mdl = max(fmdl.nodes(:,1))-min(fmdl.nodes(:,1)); % Works for circular models only
    ell_per = sum(sum( perim(2:end-1,2:end-1) )) *sz_mdl/size(orig_inhomo,1);
else
    % Knud Thomsen formula, http://planetcalc.com/149/
    p = 1.6075;
    a = ell_sz(1)^p; b = ell_sz(2)^p; c = ell_sz(3)^p;
    ell_per = 4*pi*((a*b + a*c + b*c)/3)^(1/p);
end

err_area_pct = area_err / area_frac; % pct of the semiellipse considered as error (area)

avg_err_dist = area_err / (ell_per*0.5);

end

function [slice_, chk_dir] = process_args(varargin)
slice_.do = false;
chk_dir.do = false;
if nargin==1 && isnumeric(varargin{1}) % backward compatibility
    slice_.do = true;
    slice_.slice = varargin{1};
    return;
end

k=1;
while(k<=nargin)
    if ischar(varargin{k})
        switch (varargin{k})
            case 'slice'
                slice_.do = true;
                slice_.slice = varargin{k+1};
                k = k+2;
            case 'chk_dir'
                if (k<nargin) && islogical(varargin{k+1})
                    chk_dir.do = varargin{k+1}; k = k+2;
                else
                    chk_dir.do = true; k = k+1;
                end
        end
    else % Unknown input, ignored
        k = k+1;
        warning('Input parameter not understood, ignoring');
    end
end
end
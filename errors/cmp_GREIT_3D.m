function params = cmp_GREIT_3D(imgs, good_img)
% CMP_GREIT_3D: calculate GREIT figures of merit for 3D images
% USAGE:
% params = eval_GREIT_fig_merit(imgs, xyzr_pt)
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Resolution
%  params(4,:) = Shape Deformation
%  params(5,:) = Ringing
%
%  imgs:    a sequence of eidors images of single point targets
%  good_img: The target image

% (C) 2015/10/26 Sebastien Martin

mdl = imgs(1).fwd_model;
len_x = abs( max(mdl.nodes(:,1)) - min(mdl.nodes(:,1)) );
len_z = abs( max(mdl.nodes(:,3)) - min(mdl.nodes(:,3)) );
n_slices = round( 32*(len_z/len_x) );
% Just in case one day EIDORS uses 'npz'
% [imgs(:).fwd_model.mdl_slice_mapper]= deal(struct('npx',32,'npy',32,'npz',n_slices));
minZ = min(mdl.nodes(:,3));
maxZ = max(mdl.nodes(:,3));
levels = [Inf(n_slices,2), linspace(minZ, maxZ, n_slices)'];
[imgs(:).calc_slices.levels] = levels;

% good_img.fwd_model.mdl_slice_mapper = struct('npx',32,'npy',32,'npz',n_slices);
good_img.calc_slices.levels = levels;

imgs = calc_slices(imgs);
map = ~isnan(squeeze(imgs(:,:,1,:))); %assume all imgs are the same shape
imgs(isnan(imgs)) = 0;
sz = size(imgs); sz = sz([1 2 4]);
[x,y,z,bb_min,bb_max]=prepare_grid(sz,mdl);

good_img = calc_slices(good_img);
good_img(isnan(good_img)) = 0;
[xreal,yreal,zreal, eqc_real, qmi_real,good_img] = calc_cofg( squeeze(good_img(:,:,1,:)),map,x,y,z);
real_coords = [xreal; yreal; zreal];

N_imgs = size(imgs,3);
for i= 1:N_imgs
    [xmean,ymean,zmean,equiv_circ,qmi,img] = calc_cofg(squeeze(imgs(:,:,i,:)),map,x,y,z);
    params(1,i) = calc_amplitude( img );
    params(2,i) = calc_posn_error( qmi, xmean, ymean, zmean, real_coords );
    params(3,i) = calc_diff_resolution( qmi, qmi_real, map );
    params(4,i) = calc_shape_deform( qmi, qmi_real );
    params(5,i) = calc_ringing( img, qmi );
end

% TODO: Fix this when we start to care about units
ctr = bb_min + 0.5*(bb_max-bb_min);
r = max(0.5*(bb_max-bb_min));
if N_imgs > 10 % doesn't make sense to normalize otherwise
    ctr_pts = sum((xyzr_pt(1:mdl_dim(mdl(1)),:)-repmat(ctr',1,size(xyzr_pt,2))).^2) < (0.05*r)^2;
    if any(ctr_pts)
        params(1,:) = params(1,:)/mean(params(1,ctr_pts));
    else
        eidors_msg('eval_GREIT_fig_merit: no centre points found to normalize',1);
    end
end


function ampl = calc_amplitude(img)
ampl = sum(img(:));

function pe   = calc_posn_error(qmi, xmean, ymean, zmean, real_coords)
% This definition allows + and - PE, but can also give zero in unexpected places
pe = sqrt(sum(real_coords.^2)) - sqrt( xmean^2 + ymean^2 + zmean^2 );
% This definition gives the absolute PE, but can't be negative
pe = sqrt((real_coords(1,:)-xmean).^2 + (real_coords(2,:)-ymean).^2 + (real_coords(3,:)-zmean).^2);

function dres  = calc_diff_resolution(qmi, qmi_real, map)
res = sqrt( sum(qmi(:)) / sum(map(:)));
res_real = sqrt( sum(qmi_real(:)) / sum(map(:)));
dres = abs(res-res_real);

function sd  = calc_shape_deform(qmi, qmi_real)
not_circ= qmi & ~qmi_real;
sd = sum(not_circ(:))/sum(qmi(:));

function rr = calc_ringing(img, qmi )
ring_part =  img .* ( (img<0) & ~qmi);
rr = -sum( ring_part(:) )/sum( img(:).*qmi(:) );

function [x,y,z,bb_min,bb_max]=prepare_grid(sz,mdl)
% bounding box
bnd = unique(mdl.boundary);
bb_min = min(mdl.nodes(bnd,:));
bb_max = max(mdl.nodes(bnd,:));

[x,y,z]=meshgrid(linspace(bb_min(1),bb_max(1),sz(1)),...
    linspace(bb_min(2),bb_max(2),sz(2)), linspace(bb_min(3),bb_max(3),sz(3)));


function [xmean,ymean,zmean,equiv_circ,qmi,img] = calc_cofg(img,map,x,y,z);
%  if abs(max(img(:))) < abs(min(img(:))); img= -img; end
sz = size(img,1);
qmi = calc_ha_set( img, 0.25 );
if sum(img(:) & qmi(:))<0 ;
    error('problem in CofG calculation');
end

pix_sz = (max(x(:))-min(x(:)))*(max(y(:))-min(y(:)))*(max(z(:))-min(z(:))) /numel(img);

%map = x.^2+y.^2<1.1;
qmi = qmi.*map; img = img.*map;

%  qmi = qmi .* img;  %USE THE IMAGE AMPLITUDE FOR THE CALCULATION

ss_qmi = sum(qmi(:));
xmean =  sum( sum(sum( (qmi.*x) )))/ss_qmi; % centre of gravity
ymean =  sum( sum(sum( (qmi.*y) )))/ss_qmi;
zmean =  sum( sum(sum( (qmi.*z) )))/ss_qmi;
equiv_circ = (x-xmean).^2 + (y-ymean).^2 + (z-zmean).^2 < pix_sz*ss_qmi/pi;

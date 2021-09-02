function params = eval_GREIT_3D(imgs, xyzr_pt)
% EVAL_GREIT_3D: calculate GREIT figures of merit for 3D images
% USAGE:
% params = eval_GREIT_fig_merit(imgs, xyzr_pt)
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Resolution
%  params(4,:) = Shape Deformation
%  params(5,:) = Ringing
%
%  imgs:    a sequence of eidors images of single point targets
%  xyzr_pt: [x;y;z;radius] of each images point

% (C) 2015 Sebastien Martin
% Based on eval_GREIT_fig_merit from Andy Adler

mdl = imgs.fwd_model;
[imgs(:).fwd_model.mdl_slice_mapper.npx]= deal(32);
[imgs(:).fwd_model.mdl_slice_mapper.npy]= deal(32);
[imgs(:).fwd_model.mdl_slice_mapper.npz]= deal(32); % Just in case one day EIDORS use this
minZ = min(mdl.nodes(:,3)); maxZ = max(mdl.nodes(:,3));
levels = [Inf(32,2), linspace(minZ, maxZ, 32)'];
[imgs(:).calc_slices.levels] = levels;

imgs = calc_slices(imgs);
map = ~isnan(squeeze(imgs(:,:,1,:))); %assume all imgs are the same shape
imgs(isnan(imgs)) = 0;
sz = size(imgs,1);
[x,y,z,bb_min,bb_max]=prepare_grid(sz,mdl);

N_imgs = size(imgs,3);
for i= 1:N_imgs
    [xmean,ymean,zmean,equiv_circ,qmi,img] = calc_cofg(squeeze(imgs(:,:,i,:)),map,x,y,z);
    params(1,i) = calc_amplitude( img );
    params(2,i) = calc_posn_error( qmi, xmean, ymean, zmean, xyzr_pt(1:3,i) );
    params(3,i) = calc_resolution( qmi, map );
    params(4,i) = calc_shape_deform( qmi, equiv_circ );
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

function pe   = calc_posn_error(qmi, xmean, ymean, zmean, xyz)
% This definition allows + and - PE, but can also give zero in unexpected places
pe = sqrt(sum(xyz.^2)) - sqrt( xmean^2 + ymean^2 + zmean^2 );
% This definition gives the absolute PE, but can't be negative
pe = sqrt((xyz(1,:)-xmean).^2 + (xyz(2,:)-ymean).^2 + (xyz(3,:)-zmean).^2);

function res  = calc_resolution(qmi, map)
res = sqrt( sum(qmi(:)) / sum(map(:)));

function sd  = calc_shape_deform(qmi, equiv_circ)
not_circ= qmi & ~equiv_circ;
sd = sum(not_circ(:))/sum(qmi(:));

function rr = calc_ringing(img, qmi );
ring_part =  img .* ( (img<0) & ~qmi);
rr = -sum( ring_part(:) )/sum( img(:).*qmi(:) );

function [x,y,z,bb_min,bb_max]=prepare_grid(sz,mdl)
% bounding box
bnd = unique(mdl.boundary);
bb_min = min(mdl.nodes(bnd,:));
bb_max = max(mdl.nodes(bnd,:));

[x,y,z]=meshgrid(linspace(bb_min(1),bb_max(1),sz),...
    linspace(bb_min(2),bb_max(2),sz), linspace(bb_min(3),bb_max(3),sz));


function [xmean,ymean,zmean,equiv_circ,qmi,img] = calc_cofg(img,map,x,y,z);
%  if abs(max(img(:))) < abs(min(img(:))); img= -img; end
sz = size(img,1);
%    qmi = calc_hm_set( img, 0.25 );
qmi = calc_ha_set( img, 0.25 );
if sum(img(:) & qmi(:))<0
    error('problem in CofG calculation');
end

pix_sz = (max(x(:))-min(x(:)))*(max(y(:))-min(y(:)))*(max(z(:))-min(z(:))) /numel(img);
pix_sz_tmp = (max(x(:))-min(x(:)))*(max(y(:))-min(y(:)));
num_elem_lay = squeeze( cellfun(@(x) numel(x), num2cell(img,[1 2])) );
pix_sz_2D = arrayfun(@(x) pix_sz_tmp/x, num_elem_lay);

%map = x.^2+y.^2<1.1;
qmi = qmi.*map; img = img.*map;

%  qmi = qmi .* img;  %USE THE IMAGE AMPLITUDE FOR THE CALCULATION

ss_qmi = sum(qmi(:));
ss_qmi_2D = squeeze( cellfun(@(x) sum(x(:)), num2cell(qmi,[1 2])) );

xmean =  sum( sum(sum( (qmi.*x) )))/ss_qmi; % centre of gravity
ymean =  sum( sum(sum( (qmi.*y) )))/ss_qmi;
zmean =  sum( sum(sum( (qmi.*z) )))/ss_qmi;
%    equiv_circ = (x-xmean).^2 + (y-ymean).^2 + (z-zmean).^2 < pix_sz*ss_qmi/pi;
for k=1:1:sz
    equiv_circ(:,:,k) = (x(:,:,k)-xmean).^2 + (y(:,:,k)-ymean).^2 < pix_sz_2D(k)*ss_qmi_2D(k)/pi;
end
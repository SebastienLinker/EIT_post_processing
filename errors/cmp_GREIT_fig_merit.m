function params = cmp_GREIT_fig_merit(imgs, xyzra, t_type)
% CMP_GREIT_FIG_MERIT: compare GREIT figures of merit for images
% Built-in function eval_GREIT_fig_merit compares the EIT images to a
% circular target, and therefore is not suitable for more complicated
% targets. This function actually compares the resulting images to the
% expected image, generated manually, and thus is suitable for any kind of
% target.
%
% USAGE:
% params = eval_GREIT_fig_merit(imgs, good_img)
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Difference of resolution
%  params(4,:) = Shape Deformation
%  params(5,:) = Ringing
%
%  imgs:    a sequence of eidors images of single point targets
%  good_img: The expected result

% (C) 2015/07/18 Sebastien Martin

debug = false;

if isstruct(xyzra) && isfield(xyzra,'type') && ischar(xyzra.type) && strcmp(xyzra.type,'image')
    use_depr_def = true;
    good_img = xyzra;
    eidors_msg(['You are using a deprecated call to ',mfilename,'. SD may be wrong!'],3);
else
    use_depr_def = false;
    if ~iscell(t_type) t_type = {t_type}; end
    if length(t_type)==1 && size(xyzra,2)>1
        t_type = repmat(t_type,1,size(xyzra,2));
    end
    ratio_x = max(imgs.fwd_model.mdl_slice_mapper.x_pts) - min(imgs.fwd_model.mdl_slice_mapper.x_pts);
    ratio_y = max(imgs.fwd_model.mdl_slice_mapper.y_pts) - min(imgs.fwd_model.mdl_slice_mapper.y_pts);
end

mdl = imgs.fwd_model;
imgs = calc_slices(imgs);
map = ~isnan(squeeze(imgs(:,:,1))); %assume all imgs are the same shape
imgs(isnan(imgs)) = 0;
sz = size(imgs);
[x,y,bb_min,bb_max]=prepare_grid(sz,mdl);

N_imgs = size(imgs,3);

if use_depr_def
    good_img = calc_slices(good_img);
    good_img(isnan(good_img)) = 0;
    [xreal,yreal,~, qmi_real,good_img] = calc_cog(good_img,map,x,y);
    if debug; figure; show_slices(good_img); title('expected image'); end
    xreal = repmat(xreal,1,N_imgs); yreal = repmat(yreal,1,N_imgs);
else
    fact = mean( (abs(bb_max)+abs(bb_min))/2 );
    xreal = xyzra(1,:); yreal = xyzra(2,:);
    if size(xyzra,1)>4; add_params = xyzra(5:end,:); xyzra = xyzra(1:4,:);
    else add_params = zeros(1,N_imgs);
    end
    RES_real = get_real_RES( xyzra(4,:)/fact, t_type);
end


for i= 1:N_imgs
    if use_depr_def
        [xmean,ymean,~, qmi,img] = calc_cog(imgs(:,:,i),map,x,y);
    else
        if strcmp(t_type{i},'ellipse');
            add_params(1,i) = add_params(1,i)/ratio_x; add_params(2,i) = add_params(2,i)/ratio_y;
        end
        [xmean,ymean, eq_fig, qmi,img] = calc_cog(imgs(:,:,i),map,x,y, t_type{i}, add_params(:,i));
    end
    if debug; figure; show_slices(img); title('expected image'); end
    params(1,i) = calc_amplitude( img );
    params(2,i) = calc_posn_error( qmi, xmean, ymean, xreal(i), yreal(i) );
    if use_depr_def
        params(3,i) = calc_diff_resolution( qmi, qmi_real, map );
        params(4,i) = calc_shape_deform( qmi, qmi_real );
    else
        params(3,i) = calc_diff_resolution( qmi, RES_real(i), map );
        params(4,i) = calc_shape_deform( qmi, eq_fig );
    end
    params(5,i) = calc_ringing( img, qmi );
end

% TODO: Fix this when we start to care about units
ctr = bb_min + 0.5*(bb_max-bb_min);
r = max(0.5*(bb_max-bb_min));
if N_imgs > 10 % Need many images to better normalize
    ctr_pts = sum((good_img(1:mdl_dim(mdl(1)),:)-repmat(ctr',1,size(good_img,2))).^2) < (0.05*r)^2;
    if any(ctr_pts)
        params(1,:) = params(1,:)/mean(params(1,ctr_pts));
    else
        eidors_msg('eval_GREIT_fig_merit: no centre points found to normalize',1);
    end
end


function ampl = calc_amplitude(img)
ampl = sum(img(:));

function pe   = calc_posn_error(qmi, xmean, ymean, xreal, yreal)
% This definition allows + and - PE, but can also give zero in unexpected places
pe = sqrt( xreal^2 + yreal^2 ) - sqrt( xmean^2 + ymean^2);
% This definition gives the absolute PE, but can't be negative
%  pe = sqrt((xy(1,:) - xmean).^2 + ...
%            (xy(2,:) - ymean).^2);

function dres  = calc_diff_resolution(qmi, qmi_real, map)
res = sqrt( sum(qmi(:)) / sum(map(:)));
if numel(qmi_real)==1
    dres = abs(res-qmi_real);
else % Deprecated interface
    res_real = sqrt( sum(qmi_real(:)) / sum(map(:)));
    dres = abs(res-res_real);
end

function sd  = calc_shape_deform(qmi, qmi_real)
not_circ= qmi & ~qmi_real;
sd = sum(not_circ(:))/sum(qmi(:));

function rr = calc_ringing(img, qmi );
ring_part =  img .* ( (img<0) & ~qmi);
rr = -sum( ring_part(:) )/sum( img(:).*qmi(:) );

function [x,y,bb_min,bb_max]=prepare_grid(sz,mdl)
% bounding box
bnd = unique(mdl.boundary);
bb_min = min(mdl.nodes(bnd,:));
bb_max = max(mdl.nodes(bnd,:));

[x,y]=ndgrid(linspace(bb_min(1),bb_max(1),sz(1)),linspace(bb_min(2),bb_max(2),sz(2)));


function [xmean,ymean, equiv_fig, qmi,img] = calc_cog(img,map,x,y, t_type, varargin);
qmi = calc_hm_set( img, 0.25 );
debug = evalin('caller','debug');
if sum(img(:) & qmi(:))<0
    error('problem in CofG calculation');
end

pix_sz = (max(x(:)) - min(x(:))) *( max(y(:)) - min(y(:))) /numel(img);
qmi = qmi.*map; img = img.*map;

ss_qmi = sum(qmi(:));
xmean =  sum(sum( (qmi.*x) ))/ss_qmi; % centre of gravity
ymean =  sum(sum( (qmi.*y) ))/ss_qmi;
if nargin>5
    switch(t_type)
        case {'circle','sphere','cylinder'} % This is 2D, no diff
            equiv_fig = (x-xmean).^2 + (y-ymean).^2 < pix_sz*ss_qmi/pi;
        case 'triangle'
            % X and Y arrays appear to be switched at some point before
            angle = varargin{1};
            equiv_circ = (x-xmean).^2 + (y-ymean).^2 < pix_sz*ss_qmi/pi;
            est_side = sqrt( 4*sum(equiv_circ(:))/sqrt(3) ); % Estimate length of one side?
            xm_d = ((ymean+1)/2)*size(qmi,1);
            ym_d = ((xmean+1)/2)*size(qmi,2);
            if debug
                figure; show_slices(img); axis on; hold on
                plot(xm_d,ym_d,'g*'); title('Image');
                drawTriangle(xm_d, ym_d, est_side, varargin{1});
            end
            equiv_fig = draw_triangle(xm_d, ym_d, y,x, est_side, varargin{1});
            if debug
                figure; show_slices(equiv_fig); axis on; hold on
                plot(xm_d,ym_d,'g*'); title('Shape (same area as image)');
                drawTriangle(xm_d, ym_d, est_side, varargin{1});
            end
        case 'ellipse'
            r = varargin{1};
            equiv_area = ss_qmi;
            orig_area = pi * r(1) * r(2) * size(qmi,1) * size(qmi,2);
            fact = equiv_area/orig_area;
            r = r*sqrt(fact);
            sz_mdl = max(x(:)) - min(x(:));
            equiv_fig = (x-xmean).^2/(r(2).^2*size(qmi,2)) + (y-ymean).^2/(r(1).^2*size(qmi,1)) < ...
                pix_sz*sz_mdl;
            if debug
                sz_mdl = [min(x(:)) max(x(:)); min(y(:)) max(y(:));...
                    max(x(:))-min(x(:)) max(y(:))-min(y(:))];
                Xc = ((ymean - sz_mdl(1,1) )/ sz_mdl(3,1) )*size(qmi,1);
                Yc = ((xmean - sz_mdl(2,1) )/ sz_mdl(3,2) )*size(qmi,2);
                figure; show_slices(img); axis on; hold on; title('Image');
                drawEllipse([Xc Yc r(1)*size(qmi,1) r(2)*size(qmi,2)]);
                figure; show_slices(equiv_fig); axis on; hold on; title('Shape (same area as image)');
            end
    end
else
    equiv_fig = [];
end

function triangle = draw_triangle(xmean, ymean, x,y,side,alpha)
debug = evalin('caller','debug');

h = sqrt(3)/2*side;
ang = 2*pi/3;
[theta, R] = cart2pol(2/3*h, 0);
% Define the position of the points
[pts(1,1), pts(1,2)] = pol2cart(theta+alpha, R);
[pts(2,1), pts(2,2)] = pol2cart(theta+ang+alpha, R);
[pts(3,1), pts(3,2)] = pol2cart(theta-ang+alpha, R);

sz_x = size(x,1);
pts(:,1) = (pts(:,1)+sz_x/2);% / sz_x;
sz_y = size(y,2);
pts(:,2) = (pts(:,2)+sz_y/2);% / sz_y;

pts = pts + repmat([xmean-sz_x/2,ymean-sz_y/2],3,1); %shift

if debug
    plot(pts(1,1),pts(1,2),'g*'); plot(pts(2,1),pts(2,2),'g*'); plot(pts(3,1),pts(3,2),'g*');
end

% Now go back to something between -1 and 1
xmean = 2*xmean/sz_x-1; ymean = 2*ymean/sz_y-1;
pts(:,1) = 2* pts(:,1)/sz_x -1; pts(:,2) = (2*pts(:,2)/sz_y-1);

sl(1) = (pts(2,2)-pts(1,2))./(pts(2,1)-pts(1,1)); % Slope A
sl(2) = (pts(3,2)-pts(2,2))./(pts(3,1)-pts(2,1)); % Slope B
sl(3) = (pts(1,2)-pts(3,2))./(pts(1,1)-pts(3,1)); % Slope C


% To do, rewrite this part properly
if alpha<deg2rad(-60)
    eqn{1} = (y-pts(1,2)) > (x-pts(1,1))*sl(1);
else
    eqn{1} = (y-pts(1,2)) < (x-pts(1,1))*sl(1);
end
if alpha >= 0
    eqn{2} = (y-pts(2,2)) > (x-pts(2,1))*sl(2);
else
    eqn{2} = (y-pts(2,2)) < (x-pts(2,1))*sl(2);
end
eqn{3} = (y-pts(3,2)) > (x-pts(3,1))*sl(3);
triangle = all(cat(3,eqn{:}),3);


function RES = get_real_RES(rad,t_type)
n_imgs = numel(rad);
RES = zeros(1,length(rad));
for k=1:1:n_imgs
    switch(t_type{k})
        case {'circle','sphere','cylinder'} % This is 2D, so who cares?
            area = pi*rad(k)*rad(k);
        case 'triangle'
            area = rad(k)*sqrt(3)/4;
        case 'ellipse'
            area = pi*rad(k)*rad(k);
    end
    RES(k) = sqrt( area/pi ); % Assumes circular unit mesh
end

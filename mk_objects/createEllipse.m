function [ img_data, info_obj ] = createEllipse( n_images,model,bkgnd,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,cond_r )
%CREATEELLIPSE Create ellipsoidal shape
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_size: maximal size of the object, default = 0.4
%		min_rad: for unbounded EIT, minimum radius for position
%		rad_probe: radius of the probe
%		max_pos: maximal distance from the centre (default 1)
%		min_angle: minimal angle (in degrees)
%		max_angle: maximal angle (in degrees)
%	Output
%		img_data: conductivity of the different elements in the tank, size
%		n_elements*n_images
%		info_obj: array containing size, position(X,Y) and conductivity of
%		random object, n_images*4
debug = false;
plotMatrix = false;
try
    model2 = model.fwd_model;
catch
    model2 = model;
end
model = model2;
n_elements =  size( model.elems,1 );
if (~exist('max_size','var')||isempty(max_size)) max_size = 0.4; end
if (~exist('min_rad','var')||isempty(min_rad)) min_rad = 0; end
if (~exist('rad_probe','var')||isempty(rad_probe)) rad_probe = 0; end
if (~exist('max_pos','var')||isempty(max_pos)) max_pos = 1; end
if (~exist('min_angle','var')||isempty(min_angle)) min_angle = 0; end
if (~exist('max_angle','var')||isempty(max_angle)) max_angle = 360; end
if (~exist('cond_r','var')||isempty(cond_r)) cond_r = [0 1]; end
min_angle2 = (min_angle/180)*pi; % switch to polar
max_angle2 = (max_angle/180)*pi;

%% Randomly choose the size, position and the conductivity
[inh_sizeA, inh_sizeB, inh_posX, inh_posY, inh_cond] = regenElem(n_images,max_size,min_rad,rad_probe,max_pos,min_angle2,max_angle2);
if plotMatrix
    figure; plotmatrix(inh_posX,inh_posY);
    title('Generate random ellipses within the tank: position of their centers');
end
elements = bkgnd * ones(n_elements,1);
% Added: Modify conductivity range
tot_range = sum(cond_r(:,2)-cond_r(:,1));
cond_th = [0 ((cond_r(:,2)-cond_r(:,1))/tot_range)'];
for k=2:1:size(cond_r,1)+1
    cond_th(k) = cond_th(k)+cond_th(k-1);
    tmp_elems = find( (inh_cond>=cond_th(k-1)) & (inh_cond<=cond_th(k)) );
    tmp_inhomo(tmp_elems) = inh_cond(tmp_elems)* ...
        (cond_r(k-1,2)-cond_r(k-1,1))+cond_r(k-1,1) - bkgnd; % xxx*(max-min)+min
end
inh_cond = tmp_inhomo' + bkgnd;
info_obj = [inh_sizeA inh_sizeB inh_posX inh_posY inh_cond];

%% Create elements
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    select_fcn = inline(['((x-',num2str(inh_posX(k)),')/',num2str(inh_sizeA(k)),').^2',...
        '+((y-',num2str(inh_posY(k)),')/',num2str(inh_sizeB(k)),').^2',...
        '<1'],'x','y','z');
    elements1(k,:) = bkgnd - (bkgnd-inh_cond(k)) * elem_select( img_inhomo.fwd_model, select_fcn );
    threshold = bkgnd - inh_cond(k)*0.1;
    elems_obj = find(elements1(k,:)<threshold);
    elems_no_obj = find(elements1(k,:)>=threshold);
    elements1(k, elems_obj) = inh_cond(k);
    elements1(k, elems_no_obj) = bkgnd;
    img_inhomo_rec(k) = img_inhomo;
    img_inhomo_rec(k).elem_data = elements1(k,:)';
    if debug
        figure; show_fem(img_inhomo_rec(k),[1,0,0]); title('elements1>0, elements2<1');
    end
end

%% Count the number of homogeneous images
%  We need to do it again
to_change = [];
for k=1:1:n_images
    if( min(img_inhomo_rec(k).elem_data)==max(img_inhomo_rec(k).elem_data) )
        to_change = [to_change k];
    end
end
nb_to_change = numel(to_change);
if (nb_to_change > 0)
    [tmp_new_arr, new_info] = createEllipse(nb_to_change,model,bkgnd,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,cond_r);
    info_obj(to_change,:) = new_info;
    for k=1:1:nb_to_change
        curr = to_change(k);
        img_inhomo_rec(curr).elem_data = tmp_new_arr(:,k);
    end
end

%% Create resulted matrix
img_data = zeros( n_elements, n_images );
for k=1:1:n_images
    img_data(:,k) = img_inhomo_rec(k).elem_data;
end

end

%% Regenerate new elements when we leave the circle
% n_to_regen = number of images we need to regenerate, we need to move the
% object to another random position
function [a b X Y conduc] = regenElem (n_to_regen,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,correction)
if (~exist('correction','var')) correction = true; end
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,5);
a = inhomogeneities(:,1) * (max_size-min_rad) + min_rad;
b = inhomogeneities(:,5) * (max_size-min_rad) + min_rad;
X = (inhomogeneities(:,2)*2-1)*max_pos;
Y = (inhomogeneities(:,3)*2-1)*max_pos;
conduc = inhomogeneities(:,4);

% Ensure we don't leave the circle
fix = find( sqrt(X.^2+Y.^2)>max_pos );
% Ensure the center is not where the probe is
% We might want to do this differently, why only the center? Code is
% easier to understand like that
fix = [ fix; find(sqrt(X.^2+Y.^2)<repmat(rad_probe,n_to_regen,1)) ];
fix = unique(fix);
if (numel(fix>0))
    [a(fix) b(fix) X(fix) Y(fix) conduc(fix)] = regenElem(numel(fix), max_size, min_rad, rad_probe, max_pos, min_angle, max_angle, false);
end
% Fix the angle
% Should be done ONCE (within a recursive function)
if correction
    [theta, rho] = cart2pol(X,Y);
    theta = (theta+pi)*((max_angle-min_angle)/(2*pi))+min_angle -pi;
    [X, Y] = pol2cart(theta, rho);
end
end

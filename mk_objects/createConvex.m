function [ img_data, info_obj ] = createConvex( n_images,model,bkgnd,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle, cond_r)
%CREATECONVEX Create convex elements
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_size: maximal size of the object, default = 0.4
%		min_rad: for unbounded EIT, minimal distance from center
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
[inhomo_size, inhomo_positionX, inhomo_positionY, inhomo_conduct] = regenElem(n_images,max_size,min_rad,rad_probe,max_pos,min_angle2,max_angle2);
if plotMatrix
    figure; plotmatrix(inhomo_positionX,inhomo_positionY);
    title('Generate random objects within the tank: position of their centers');
end
elements = bkgnd * ones(n_elements,1);
if (numel(unique(cond_r(:)))==1)
    inhomo_conduct = inhomo_conduct*(cond_r(1)-cond_r(1))+cond_r(1);
else
    tot_range = sum(cond_r(:,2)-cond_r(:,1));
    cond_th = [0 ((cond_r(:,2)-cond_r(:,1))/tot_range)'];
    for k=2:1:size(cond_r,1)+1
        cond_th(k) = cond_th(k)+cond_th(k-1);
        tmp_elems = find( (inhomo_conduct>=cond_th(k-1)) & (inhomo_conduct<=cond_th(k)) );
        tmp_inhomo(tmp_elems) = inhomo_conduct(tmp_elems)* ...
            (cond_r(k-1,2)-cond_r(k-1,1))+cond_r(k-1,1) - bkgnd; % xxx*(max-min)+min
    end
    inhomo_conduct = tmp_inhomo' + bkgnd;
end

info_obj = [inhomo_size inhomo_positionX inhomo_positionY inhomo_conduct];

%% Make the object convex
% Generate 2 images
% Thresholding algorithm
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    select_fcn = inline(['(x-',num2str(inhomo_positionX(k)),').^2+(y-', ...
        num2str(inhomo_positionY(k)),').^2<',num2str(inhomo_size(k)),'^2'],'x','y','z');
    
    elements1(k,:) = bkgnd - inhomo_conduct(k) * elem_select( img_inhomo.fwd_model, select_fcn );
    threshold = bkgnd - inhomo_conduct(k)*0.1;
    % Change 2014/07/08: variable conductivity range
    elems_obj = find(elements1(k,:)<threshold);
    elems_no_obj = find(elements1(k,:)>=threshold);
    % Important: This line have been modified, but not tested yet for
    % previous cases, backward compatibility may result, if so, please
    % use a former revision from SVN server or use the commented output
    elements1(k, elems_obj) = inhomo_conduct(k); % bkgnd - inhomo_conduct(k);
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
    [tmp_new_arr, new_info] = createConvex(nb_to_change,model,bkgnd,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,cond_r);
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
function [size,X,Y,conduc] = regenElem (n_to_regen,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,correction)
if (~exist('correction','var')) correction = true; end
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,5);
size = inhomogeneities(:,1) * (max_size-min_rad) + min_rad;
X = (inhomogeneities(:,2)*2-1)*max_pos;
Y = (inhomogeneities(:,3)*2-1)*max_pos;
conduc = inhomogeneities(:,4);
% Ensure we don't leave the circle
% Ensure the center is not where the probe is
% We might want to do this differently, why only the center? Code is
% easier to understand like that
fix = find( sqrt(X.^2+Y.^2)>max_pos );
fix = [fix; find(sqrt(X.^2+Y.^2)<rad_probe) ];
fix = unique(fix);
% Fix the angle
% Should be done ONCE (within a recursive function)
if correction
    [theta, rho] = cart2pol(X,Y);
    theta = (theta)*((max_angle-min_angle)/(2*pi))+min_angle;
    rho = rho*(max_pos-rad_probe)+rad_probe;
    [X, Y] = pol2cart(theta, rho);
end
end

function [ img_data, info_obj ] = createConvex3Dbounded( n_images,model,bkgnd,min_size,max_size,max_pos,min_angle,max_angle,Z_pos)
%CREATECONVEX3DBOUNDED Create convex elements for 3D bounded domain
%	For training purpose, all of the elements can have their centre located from 0 to 90 degrees
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		min_size: minimal size of the object, default = 0
%		max_size: maximal size of the object, default = 0.3
%		max_pos: maximal distance from the centre (default 1)
%		min_angle: minimum angle (in degrees) (default 0)
%		max_angle: maximum angle in degrees (default 360)
%		Z_pos
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
if ~exist('min_size','var') min_size = 0; end
if ~exist('max_size','var') max_size = 0.3; end
if ~exist('min_rad','var') min_rad = 0; end
if ~exist('rad_probe','var') rad_probe = 0; end
if ~exist('max_pos','var') max_pos = 1; end
if ~exist('Z_pos','var') Z_pos = [-4.3; 0]; end

%% Randomly choose the size, position and the conductivity
[size_x, size_y, size_z, inhomo_positionX, inhomo_positionY, inhomo_positionZ, inhomo_conduct] = ...
    regenElem(n_images,min_size,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,Z_pos);
if plotMatrix
    figure; plotmatrix(inhomo_positionX,inhomo_positionY);
    title('Generate random objects within the tank: position of their centers');
end
elements = bkgnd * ones(n_elements,1);
info_obj = [size_x size_y size_z inhomo_positionX inhomo_positionY inhomo_positionZ inhomo_conduct];

%% Make the object convex
% Thresholding algorithm
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    select_fcn = inline(['((x-',num2str(inhomo_positionX(k)),')/',num2str(size_x(k)),').^2+', ...
        '((y-',num2str(inhomo_positionY(k)),')/',num2str(size_y(k)),').^2+', ...
        '((z-',num2str(inhomo_positionZ(k)),')/',num2str(size_z(k)),').^2<1'],'x','y','z');
    elements1(k,:) = bkgnd - inhomo_conduct(k) * elem_select( img_inhomo.fwd_model, select_fcn );
    threshold = bkgnd - inhomo_conduct(k)*0.1;
    elements1(k, find(elements1(k,:)<threshold)) = bkgnd - inhomo_conduct(k);
    elements1(k, find(elements1(k,:)>threshold)) = bkgnd;
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
    [tmp_new_arr, new_info]= createConvex3Dbounded(nb_to_change,model,bkgnd,min_size,max_size,max_pos,min_angle,max_angle,Z_pos);
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
function [size_x, size_y, size_z, X,Y,Z, conduc] = regenElem ( n_to_regen,min_size,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,Z_pos)
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,5);
size =  inhomogeneities(:,1) * (max_size-min_size) + min_size;
size_x = size; %6*ones(n_to_regen,1)/2;
size_y = size; %9*ones(n_to_regen,1)/2;
size_z = size; %4*ones(n_to_regen,1)/2;
X = (inhomogeneities(:,2)*2-1); %.*(max_pos+size_x/2);
Y = (inhomogeneities(:,3)*2-1); %.*(max_pos+size_y/2);
Z = ( inhomogeneities(:,4).* (Z_pos(1)-Z_pos(2)) )+Z_pos(2);
conduc = inhomogeneities(:,5);

% Ensure we don't leave the circle
for k=1:1:n_to_regen
    distance(k) = sqrt( DistanceEllipse( [0,0], [X(k),Y(k),size_x(k),size_y(k),0]') );
end
dist_cent = sqrt(X.^2+Y.^2); % Use distance center to center, instead of
% minimal distance ellipse-center (should we use maximal distance?)

fix = find(dist_cent>max_pos);
% Ensure the center is not where the probe is
% We might want to do this differently, why only the center? Code is
% easier to understand like that
fix = [ fix, find(distance<rad_probe) ];
fix = unique(fix);
if (numel(fix>0))
    [size_x(fix) size_y(fix) size_z(fix) X(fix) Y(fix) Z(fix) conduc(fix)] = regenElem(numel(fix),min_size,max_size,min_rad,rad_probe,max_pos,min_angle,max_angle,Z_pos);
end

end

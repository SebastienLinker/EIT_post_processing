function [ img_data, info_obj ] = createEllipse3D( n_images,model,bkgnd,max_size,min_rad,rad_probe,max_pos)
%CREATECONVEX Create convex elements
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_size: maximal size of the object, default = 0.4
%		min_rad: for unbounded EIT, radius of the probe
%	`	rad_probe: radius of the probe
%		max_pos: maximal distance from the centre (default 1)
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
if ~exist('max_size','var') max_size = 0.4; end
if ~exist('min_rad','var') min_rad = 0; end
if ~exist('rad_probe','var') rad_probe = 0; end
if ~exist('max_pos','var') max_pos = 1; end

%% Randomly choose the size, position and the conductivity
[inhomoA, inhomoB, inhomoC, inhomo_positionX, inhomo_positionY, inhomo_conduct] = regenElem(n_images,max_size,min_rad,rad_probe,max_pos);
if plotMatrix
    figure; plotmatrix(inhomo_positionX,inhomo_positionY);
    title('Generate random objects within the tank: position of their centers');
end
elements = bkgnd * ones(n_elements,1);
info_obj = [inhomoA inhomoB inhomoC inhomo_positionX inhomo_positionY inhomo_conduct];

%% Make the object convex
% Generate 2 images
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    select_fcn = inline(['((x-',num2str(inhomo_positionX(k)),')/',num2str(inhomoA(k)),').^2',...
        '+((y-',num2str(inhomo_positionY(k)),')/',num2str(inhomoB(k)),').^2',...
        '+((z+2.4)/',num2str(inhomoC(k)),').^2',...
        '<1'],'x','y','z');
    elements1(k,:) = bkgnd - inhomo_conduct(k) * elem_select( img_inhomo.fwd_model, select_fcn );
end

% Adjust the shape
orig_elem = zeros(n_elements, n_images); new_elem = zeros(n_elements, n_images);
for k = 1:1:n_images
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
    [tmp_new_arr, new_info] = createEllipse3D(nb_to_change,model,bkgnd,max_size,min_rad,rad_probe,max_pos);
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
function [a b c X Y conduc] = regenElem ( n_to_regen, max_size, min_rad, rad_probe, max_pos)
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,6);
size = inhomogeneities(:,1) * (max_size-min_rad) + min_rad;
a = 0.4; b = 0.6; c = 0.9;
X = (inhomogeneities(:,2)*2-1);
Y = (inhomogeneities(:,3)*2-1);
conduc = 0.5; %inhomogeneities(:,4);
% Ensure we don't leave the circle
fix = find(sqrt(X.^2+Y.^2)+size>1);
fix = [fix; find( sqrt(X.^2+Y.^2)>max_pos )];
% Ensure the center is not where the probe is
% We might want to do this differently, why only the center? Code is
% easier to understand like that
fix = [ fix; find(sqrt(X.^2+Y.^2)<repmat(rad_probe,n_to_regen,1)) ];
fix = unique(fix);
if (numel(fix>0))
    [a(fix) b(fix) c(fix) X(fix) Y(fix) conduc(fix)] = regenElem(numel(fix), max_size, min_rad, rad_probe, max_pos);
    r = 0.15;
    X = a/2+r+rad_probe; Y=0;
end
end

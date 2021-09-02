function [ img_data, info_obj, eqn ] = createTriangle( n_images,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r )
%CREATETRIANGLE Create an equilateral triangle
% Inspired from function createEllipseGeneral
%	Important: The input is the size of the inner circle of the triangle
%	while the output size is the length of the side of the triangle
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_sz: maximal size of the inner circle, default: 0.4
%		min_sz: minimum size, default: 0
%		rad_probe: radius of the probe
%		max_pos: maximal distance from the centre (default 1)
%				3D unbounded: We may need s pecific Z-axis: In this case,
%				please use a cell array (2*1) or a 3*1 vector with:
%				Maximal distance (default 1)
%				Position of Z-axis (minimum and maximum)
%				For fixed Z-axis, a 2*1 vector is also ok
%		rot_angle: The rotation angle (2*1 array with minimum and maximum, in degrees, default [-90 0])
%		pos_angle: Position angle, (2*1 array with minimum and maximum, in degrees, default [0 360])
%		cond_r:
%	Output
%		img_data: conductivity of the different elements in the tank, size
%		n_elements*n_images
%		info_obj: array containing size of one side, position(X,Y,Z), rotation angles, and conductivity of
%		random object, n_images*6 (2D) or n_images*9 (3D)
%		eqn: The set of equations of the rectangle
%
%	(C) 2015/05/20 Sebastien Martin
%

debug = false;
plotMatrix = false;
try
    model2 = model.fwd_model;
catch
    model2 = model;
end
model = model2;
is2D = (size(model.nodes,2)==2);
n_elements =  size( model.elems,1 );
if (~exist('max_sz','var')||isempty(max_sz)); max_sz = 0.4; end
if (~exist('min_sz','var')||isempty(min_sz)); min_sz = 0; end
if (~exist('rad_probe','var')||isempty(rad_probe)); rad_probe = 0; end
if (~exist('max_pos','var')||isempty(max_pos)); max_pos = 1; end
if iscell(max_pos); max_pos=cell2mat(max_pos); end
if length(max_pos)==2; max_pos=[max_pos max_pos(2)]; end
if (~exist('rot_angle','var')||isempty(rot_angle)); rot_angle = [-90 0]; end % We can set this to -120 0
if (~exist('pos_angle','var')||isempty(pos_angle)); pos_angle = [0 360]; end
if (~exist('cond_r','var')||isempty(cond_r)); cond_r = [0 1]; end

rot_angle( rot_angle>0 ) = -rot_angle(rot_angle>0);
if (length(max_sz)==1); max_sz=[max_sz, max_sz]; end
if (length(min_sz)==1); min_sz=[min_sz, min_sz]; end
if (length(rot_angle)==1); rot_angle=[rot_angle, rot_angle]; end
if abs(diff(rot_angle))>120; rot_angle = [-120 0]; end
if length(pos_angle)==1; pos_angle=[pos_angle, pos_angle]; end
min_ang = (pos_angle(1)/180)*pi; % switch to polar
max_ang = (pos_angle(2)/180)*pi;


%% Randomly choose the size, position and the conductivity
if is2D
    [tri_side, rot1, pX,pY, inh_cond] = regenElem2D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
else
    [tri_side, rot1,rot2, pX,pY,pZ, inh_cond] = regenElem3D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
end
if plotMatrix
    figure; plotmatrix(pX,pY);
    title('Generate random parallelograms within the tank: position of their CoGs');
end
elements = bkgnd * ones(n_elements,1);
% Modify conductivity range
if (numel(unique(cond_r(:)))==1)
    inh_cond = cond_r(1); % xxx*(max-min)+min
else
    tot_range = sum(cond_r(:,2)-cond_r(:,1));
    cond_th = [0 ((cond_r(:,2)-cond_r(:,1))/tot_range)'];
    for k=2:1:size(cond_r,1)+1
        cond_th(k) = cond_th(k)+cond_th(k-1);
        tmp_elems = find( (inh_cond>=cond_th(k-1)) & (inh_cond<=cond_th(k)) );
        tmp_inhomo(tmp_elems) = inh_cond(tmp_elems)* ...
            (cond_r(k-1,2)-cond_r(k-1,1))+cond_r(k-1,1) - bkgnd; % xxx*(max-min)+min
    end
    inh_cond = tmp_inhomo' + bkgnd;
end
% inhomo_conduct = inhomo_conduct*(cond_r(2)-cond_r(1))+cond_r(1)-bkgnd; % xxx*(max-min)+min
if is2D
    info_obj = [tri_side rad2deg(rot1) pX pY inh_cond zeros(n_images,1)];
else
    info_obj = [tri_side rad2deg(rot1) rad2deg(rot2) pX pY pZ inh_cond];
end

%% Create elements
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    % x'=x*cos(theta)+y*sin(theta), y'=y*cos(theta)-x*sin(theta)
    if is2D
        eqn = eqn_2D(pX(k),pY(k), rot1(k), tri_side(k));
    else
        eqn = eqn_3D(pX(k),pY(k),pZ(k), rot1(k),rot2(k), tri_side(k));
    end
    elements1(k,:) = bkgnd - (bkgnd-inh_cond(k)) * elem_select( img_inhomo.fwd_model, eqn );
    if inh_cond<bkgnd
        threshold = bkgnd - (bkgnd-inh_cond(k))*0.1;
        elems_obj = find(elements1(k,:)<threshold);
        elems_no_obj = find(elements1(k,:)>=threshold);
    else
        threshold = bkgnd - (bkgnd-inh_cond(k))*0.1;
        elems_obj = find(elements1(k,:)>threshold);
        elems_no_obj = find(elements1(k,:)<=threshold);
    end
    elements1(k, elems_obj) = inh_cond(k); % bkgnd - inhomo_conduct(k);
    elements1(k, elems_no_obj) = bkgnd;
    img_inhomo_rec(k) = img_inhomo;
    img_inhomo_rec(k).elem_data = elements1(k,:)';
    if debug
        figure; show_fem(img_inhomo_rec(k),[1,0,0]); title('elements1>0, elements2<1');
        hold on; drawTriangle( pX(k), pY(k), tri_side(k), rot1(k));
        disp(['CoG: (',num2str(pX(k)),',',num2str(pY(k)),'), size: ',num2str(tri_side(k)), ...
            ' angle: ',num2str(rot1(k))]);
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
    [tmp_new_arr, new_info] = createTriangle(nb_to_change,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r);
    info_obj(to_change,:) = new_info;
    for k=1:1:nb_to_change
        curr = to_change(k);
        img_inhomo_rec(curr).elem_data = tmp_new_arr(:,k);
    end
end

%% Create resulting matrix
img_data = zeros( n_elements, n_images );
for k=1:1:n_images
    img_data(:,k) = img_inhomo_rec(k).elem_data;
end

end

%% Equations
function eqn = eqn_2D(pX, pY, alpha, tri_side)
h = sqrt(3)/2*tri_side;
ang = 2*pi/3;
[theta, R] = cart2pol(2/3*h, 0);
% Define the position of the points
[pts(1,1), pts(1,2)] = pol2cart(theta+alpha, R);
[pts(2,1), pts(2,2)] = pol2cart(theta+ang+alpha, R);
[pts(3,1), pts(3,2)] = pol2cart(theta-ang+alpha, R);

pts = pts + repmat([pX,pY],3,1); %shift
sl(1) = (pts(2,2)-pts(1,2))./(pts(2,1)-pts(1,1)); % Slope A
sl(2) = (pts(3,2)-pts(2,2))./(pts(3,1)-pts(2,1)); % Slope B
sl(3) = (pts(1,2)-pts(3,2))./(pts(1,1)-pts(3,1)); % Slope C

% To do, rewrite this part properly
if alpha<deg2rad(-60)
    eqn{1} = inline(['(y-',num2str(pts(1,2)),')>((x-',num2str(pts(1,1)),')*',num2str(sl(1)),')'],'x','y','z');
else
    eqn{1} = inline(['(y-',num2str(pts(1,2)),')<((x-',num2str(pts(1,1)),')*',num2str(sl(1)),')'],'x','y','z');
end
if alpha >= 0
    eqn{2} = inline(['(y-',num2str(pts(2,2)),')>((x-',num2str(pts(2,1)),')*',num2str(sl(2)),')'],'x','y','z');
else
    eqn{2} = inline(['(y-',num2str(pts(2,2)),')<((x-',num2str(pts(2,1)),')*',num2str(sl(2)),')'],'x','y','z');
end
eqn{3} = inline(['(y-',num2str(pts(3,2)),')>((x-',num2str(pts(3,1)),')*',num2str(sl(3)),')'],'x','y','z');
end

function eqn = eqn_3D(pX, pY, pZ, rot1, rot2, side)
eqn = eqn_2D(pX, pY, rot1, side);
end

%% Regenerate new elements
% n_to_regen = number of images we need to regenerate
function [side,rot_par,X,Y,conduc] = regenElem2D(N,max_sz,min_sz,rad_pb,max_pos,rot_ang,min_ang,max_ang,correction)
if (~exist('correction','var')); correction = true; end
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(N,5);
X = (inhomogeneities(:,1)*2-1)*max_pos;
Y = (inhomogeneities(:,2)*2-1)*max_pos;
inner_rad = inhomogeneities(:,3) * (max_sz(1)-min_sz(1)) + min_sz(1);
side = inner_rad*6/sqrt(3); % http://www.ajdesigner.com/phptriangle/equilateral_triangle_inscribed_circle_radius_r.php
rot_par = inhomogeneities(:,4) * (rot_ang(2)-rot_ang(1)) + rot_ang(1);
conduc = inhomogeneities(:,5);
% Ensure we don't leave the circle or the center is not where the probe is
% We might want to do this differently, why only the center? Code is
% easier to understand like that
% Fix the angle
[theta, rho] = cart2pol(X,Y);
theta = (theta+pi)*((max_ang-min_ang)/(2*pi))+min_ang;
rho = rho*(max_pos-rad_pb)+rad_pb;
[X, Y] = pol2cart(theta, rho);
rot_par = deg2rad(rot_par); %Convert to radians, only once
end

function [side,rot1,rot2,X,Y,Z,conduc] = regenElem3D(N,max_sz,min_sz,rad_pb,max_pos,rot_ang,min_ang,max_ang)
% Consider a straight triangle
% No rotation along the Z-axis
rot2 = zeros(N,1);
Z = zeros(N,1);
[side, rot1, X, Y, conduc] = regenElem2D(N,max_sz,min_sz,rad_pb,max_pos,rot_ang,min_ang,max_ang);
end


function [ img_data, info_obj, eqn ] = createParallelogram( n_images,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r )
%CREATEPARALLELOGRAM Create a parallelogram
% Inspired from function createEllipseGeneral
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_sz: maximal size of the object
% 			can be a 2*1 or 3*1 array with different sizes for 2 or 3 segments, default: 0.4
%		min_sz: minimum size
%			can be a 2*1 or 3*1 array with different sizes for 2 or 3 radii, default: 0
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
%		info_obj: array containing size, position(X,Y), rotation angles, and conductivity of
%		random object, n_images*6 (2D) or n_images*9 (3D)
%		eqn: The set of equations of the rectangle

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
if (~exist('rot_angle','var')||isempty(rot_angle)); rot_angle = [-90 0]; end
if (~exist('pos_angle','var')||isempty(pos_angle)); pos_angle = [0 360]; end
if (~exist('cond_r','var')||isempty(cond_r)); cond_r = [0 1]; end

if (length(max_sz)==1); max_sz=[max_sz, max_sz]; end
if (length(min_sz)==1); min_sz=[min_sz, min_sz]; end
if length(pos_angle)==1; pos_angle=[pos_angle, pos_angle]; end
min_ang = (pos_angle(1)/180)*pi; % switch to polar
max_ang = (pos_angle(2)/180)*pi;
rot_angle( rot_angle>0 ) = -rot_angle(rot_angle>0);

%% Randomly choose the size, position and the conductivity
if is2D
    [szA,szB, rot1, pX,pY, inh_cond] = regenElem2D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
else
    [szA,szB,inh_szC, rot1,rot2, pX,pY,pZ, inh_cond] = regenElem3D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
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
    info_obj = [szA szB rad2deg(rot1) pX pY inh_cond];
else
    info_obj = [szA szB inh_szC rad2deg(rot1) rad2deg(rot2) pX pY pZ inh_cond];
end

%% Create elements
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    % x'=x*cos(theta)+y*sin(theta), y'=y*cos(theta)-x*sin(theta)
    if is2D
        eqn = eqn_2D(pX(k),pY(k), rot1(k), szA(k),szB(k));
    else
        eqn = eqn_3D(pX(k),pY(k),pZ(k), rot1(k),rot2(k), szA(k),szB(k),inh_szC(k));
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
    [tmp_new_arr, new_info] = createParallelogram(nb_to_change,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r);
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
function eqn = eqn_2D(pX, pY, alpha, szA, szB)
pts = [	(szA/2) * (cos(alpha)+sin(alpha)), (szB/2) * (-cos(alpha)+sin(alpha));...
    (szA/2) * (cos(alpha)-sin(alpha)), (szB/2) * (cos(alpha)+sin(alpha)) ];
pts = [ pts; -pts]; % The opposite points
pts = pts + repmat([pX,pY],4,1); %shift
slA = (pts(2,2)-pts(1,2))/(pts(2,1)-pts(1,1)); % Slope A
slB = -1/slA; % Slope B

eqn{1} = inline(['(y-',num2str(pts(1,2)),')>((x-',num2str(pts(1,1)),')*',num2str(slA),')'],'x','y','z');
eqn{2} = inline(['(y-',num2str(pts(2,2)),')<((x-',num2str(pts(2,1)),')*',num2str(slB),')'],'x','y','z');
eqn{3} = inline(['(y-',num2str(pts(3,2)),')<((x-',num2str(pts(3,1)),')*',num2str(slA),')'],'x','y','z');
eqn{4} = inline(['(y-',num2str(pts(4,2)),')>((x-',num2str(pts(4,1)),')*',num2str(slB),')'],'x','y','z');
end

function eqn = eqn_3D(pX, pY, pZ, rot1, rot2, szA, szB, szC)
error('Unwritten code');
end

%% Regenerate new elements
% n_to_regen = number of images we need to regenerate
function [a,b,rot_par,X,Y,conduc] = regenElem2D(N,max_sz,min_sz,rad_pb,max_pos,rot_ang,min_ang,max_ang,correction)
if (~exist('correction','var')); correction = true; end
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(N,6);
X = (inhomogeneities(:,1)*2-1)*max_pos;
Y = (inhomogeneities(:,2)*2-1)*max_pos;
a = inhomogeneities(:,3) * (max_sz(1)-min_sz(1)) + min_sz(1);
b = inhomogeneities(:,4) * (max_sz(2)-min_sz(2)) + min_sz(2);
% b=a; % Square (or rhombus)
rot_par = inhomogeneities(:,5) * (rot_ang(2)-rot_ang(1)) + rot_ang(1); % Parallelogram
conduc = inhomogeneities(:,6);
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

function [a,b,c,rot1,rot2,X,Y,Z,conduc] = regenElem3D (n_to_regen,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_angle,max_angle)
error('This code hasn''t been written yet.');
% Randomly choose the size, position and the conductivity
end


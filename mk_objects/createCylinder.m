function [ img_data, info_obj, eqn ] = createCylinder( n_images,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r )
%CREATECYLINDER Creates ellipsoidal cylinder shape
% Inspired from function createEllipseGeneral
% The cylinders go along the z-axis
%   Input
%		n_images: number of images to generate
%		model: the EIDORS model
%		bkgnd: background value (homogenous tank)
%		max_sz: maximal size of the object
% 			a 2*1 array with different sizes for 2 radii, default: 0.4
%		min_sz: minimum size
%			 a 2*1 array with different sizes for 2 radii, default: 0
%		rad_probe: radius of the probe
%		max_pos: maximal distance from the centre (default 1)
%		rot_angle: The rotation angle (2*1 array with minimum and maximum, in degrees, default [0 180])
%		pos_angle: Position angle, (2*1 array with minimum and maximum, in degrees, default [0 360])
%	Output
%		img_data: conductivity of the different elements in the tank, size
%		n_elements*n_images
%		info_obj: array containing size, position(X,Y) and conductivity of
%		random object, n_images*5 or n_images*7
%		eqn: The equation of the cylinder
%
%	(C) 2016/06/19 Sebastien Martin

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
if length(pos_angle)==1; pos_angle=[pos_angle pos_angle]; end
if (~exist('rot_angle','var')||isempty(rot_angle)); rot_angle = [0 180]; end
if (~exist('pos_angle','var')||isempty(pos_angle)); pos_angle = [0 360]; end
if (~exist('cond_r','var')||isempty(cond_r)); cond_r = [0 1]; end

if (length(max_sz)==1); max_sz=[max_sz, max_sz]; end
if (length(min_sz)==1); min_sz=[min_sz, min_sz]; end
min_ang = (pos_angle(1)/180)*pi; % switch to polar
max_ang = (pos_angle(2)/180)*pi;

%% Randomly choose the size, position and the conductivity
if is2D
    [inh_szA,inh_szB, rot1, inh_pX,inh_pY, inh_cond] = regenElem2D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
else
    [inh_szA,inh_szB, rot1,rot2, inh_pX,inh_pY, inh_cond] = regenElem3D(n_images,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_ang,max_ang);
end
if plotMatrix
    figure; plotmatrix(inh_pX,inh_pY);
    title('Generate random ellipses within the tank: position of their centers');
end
elements = bkgnd * ones(n_elements,1);
% Modify conductivity range
if (numel(unique(cond_r(:)))==1)
    inh_cond = inh_cond*(cond_r(1)-cond_r(1))+cond_r(1); % xxx*(max-min)+min
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
if is2D
    info_obj = [inh_szA inh_szB rad2deg(rot1) inh_pX inh_pY inh_cond];
else
    info_obj = [inh_szA inh_szB rad2deg(rot1) rad2deg(rot2) inh_pX inh_pY inh_cond];
end

%% Create elements
img_homo = eidors_obj( 'image', 'homogeneous image', 'elem_data', elements, 'fwd_model', model);
img_inhomo = mk_image( model, elements );
elements1 = zeros(n_images, n_elements);
for k = 1:1:n_images
    if is2D
        eqn = eqn_2D(inh_pX(k),inh_pY(k), rot1(k), inh_szA(k),inh_szB(k));
    else
        eqn = eqn_3D(inh_pX(k),inh_pY(k), rot1(k),rot2(k), inh_szA(k),inh_szB(k));
    end
    elements1(k,:) = bkgnd - (bkgnd-inh_cond(k)) * elem_select( img_inhomo.fwd_model, eqn );
    if inh_cond<bkgnd
        threshold = bkgnd - inh_cond(k)*0.5;
        elems_obj = find(elements1(k,:)<threshold);
        elems_no_obj = find(elements1(k,:)>=threshold);
    else
        threshold = bkgnd + inh_cond(k)*0.1;
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
    [tmp_new_arr, new_info] = createCylinder(nb_to_change,model,bkgnd,max_sz,min_sz,rad_probe,max_pos,rot_angle,pos_angle,cond_r);
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
function eqn = eqn_2D(pX, pY, rot, szA, szB)
% x'=x*cos(theta)+y*sin(theta), y'=y*cos(theta)-x*sin(theta)
eqn = inline(['((x-',num2str(pX),')*',num2str(cos(rot)),'+',...
    '(y-',num2str(pY),')*',num2str(sin(rot)),...
    ').^2/(',num2str(szA),').^2',...
    '+((y-',num2str(pY),')*',num2str(cos(rot)),'-',...
    '(x-',num2str(pX),')*',num2str(sin(rot)),...
    ').^2/(',num2str(szB),').^2',...
    '<1'],'x','y','z');
end

function eqn = eqn_3D(pX, pY, rot1, rot2, szA, szB)
% x'=x*cos(theta)+y*sin(theta), y'=y*cos(theta)-x*sin(theta)
eqn = inline(['((x-',num2str(pX),')*',num2str(cos(rot1)),'+',...
    '(y-',num2str(pY),')*',num2str(sin(rot1)),...
    ').^2/(',num2str(szA),').^2',...
    '+((y-',num2str(pY),')*',num2str(cos(rot1)),'-',...
    '(x-',num2str(pX),')*',num2str(sin(rot1)),...
    ').^2/(',num2str(szB),').^2',...
    '<1'],'x','y','z');
end

%% Regenerate new elements
% n_to_regen = number of images we need to regenerate
function [a,b,rot,X,Y,conduc] = regenElem2D (n_to_regen,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_angle,max_angle,correction)
if (~exist('correction','var')); correction = true; end
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,6);
a = inhomogeneities(:,1) * (max_sz(1)-min_sz(1)) + min_sz(1);
b = inhomogeneities(:,5) * (max_sz(2)-min_sz(2)) + min_sz(2);
rot = inhomogeneities(:,6) * (rot_angle(2)-rot_angle(1)) + rot_angle(1);
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
% Fix the angle
if correction
    [theta, rho] = cart2pol(X,Y);
    theta = (theta+pi)*((max_angle-min_angle)/(2*pi))+min_angle;
    rho = rho*(max_pos-rad_probe)+rad_probe;
    [X, Y] = pol2cart(theta, rho);
    rot = deg2rad(rot); %Convert to radians, only once
end
end

function [a,b,rot1,rot2,X,Y,conduc] = regenElem3D (n_to_regen,max_sz,min_sz,rad_probe,max_pos,rot_angle,min_angle,max_angle)
% Randomly choose the size, position and the conductivity
inhomogeneities = rand(n_to_regen,7);
a = inhomogeneities(:,1) * (max_sz(1)-min_sz(1)) + min_sz(1);
b = inhomogeneities(:,2) * (max_sz(2)-min_sz(2)) + min_sz(2);
rot1 = inhomogeneities(:,3) * (rot_angle(2)-rot_angle(1)) + rot_angle(1);
rot2 = inhomogeneities(:,4) * (rot_angle(2)-rot_angle(1)) + rot_angle(1);
X = (inhomogeneities(:,5)*2-1);
Y = (inhomogeneities(:,6)*2-1);
conduc = inhomogeneities(:,7);

% Fix the angle
[theta, rho] = cart2pol(X,Y);
theta = (theta+pi)*((max_angle-min_angle)/(2*pi))+min_angle;
rho = rho*(max_pos(1)-rad_probe)+rad_probe;
[X, Y] = pol2cart(theta, rho);
rot1 = deg2rad(rot1); rot2 = deg2rad(rot2); %Convert to radians
end


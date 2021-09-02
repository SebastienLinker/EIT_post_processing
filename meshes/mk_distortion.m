function [ mdls_dist ] = mk_distortion( mdl, n_imgs, dist_type, varargin )
%MK_DISTORTION Make distorted meshes
%   Starting from one mesh, we can make several distortions in order to
%   solve the problem (typically the forward problem) with a distorted mesh
%
%	INPUTS:
%		mdl: Forward model
%		dist_type: a string representing the type of the distortion.
%			Possible values are:
%			-Ellipse: get an elliptical shape
%			-Complex: Complex conformal deformation
%			-Dual-Complex: 2 Complex conformal deformations
%			-Fourier: Using Fourier coefficients
%		varargin: Additional arguments, depends on dist_type
%
%	(C) 2015/06/30 Sebastien Martin


debug = false;
if (isfield(mdl,'type') && strcmp(mdl.type,'inv_model')); mdl=mdl.fwd_model; end
if ~(isfield(mdl,'type') && strcmp(mdl.type,'fwd_model')); error('Wrong parameter ''mdl'''); end
if ~ischar(dist_type); error('Please check parameter ''dist_type'''); end

params = parse_varargin( dist_type, n_imgs, varargin{:} );

switch params.type
    case 'ellipse'
        mdl =  repmat(mdl,n_imgs,1);
        mdls_dist = arrayfun(@(x,y,z) shape_deformation(x,y,z), mdl, params.X, params.Y);
    case 'conformal'
        mdl =  repmat(mdl,n_imgs,1);
        mdls_dist = arrayfun(@(x,y,z) conformal_mapping(x,y,z), mdl, params.real, params.imag);
    case 'dual_conformal'
        mdl =  repmat(mdl,n_imgs,1);
        mdls_dist = arrayfun(@(x,y,z) dual_conformal_mapping(x,y,z), mdl, params.i{1}, params.i{2});
    case 'Fourier'
        mdls_dist =  repmat(mdl,n_imgs,1);
        for k = 1:1:n_imgs
            coefs = params.fourier_coeffs;
            coefs(:,1) = coefs(:,1) .* params.fourier_max_dist_X(:,k);
            coefs(:,2) = coefs(:,2) .* params.fourier_max_dist_Y(:,k);
            mdls_dist(k) = distort_fourier(mdl,coefs);
        end
    otherwise
        error('Unknown deformation type, cannot distort the FE model');
end

if debug
    for k=1:1:n_imgs
        figure; show_fem(mdls_dist(k)); title('Distorted model');
    end
end

end

function params = parse_varargin(dist_type, n_imgs, varargin)
switch dist_type
    case 'ellipse'
        params.type = 'ellipse';
        new_x = varargin{1};
        new_y = varargin{2};
        params.X = chk_in_param(n_imgs,new_x);
        params.Y = chk_in_param(n_imgs,new_y);
    case {'complex','conformal'}
        params.type = 'conformal';
        params.real = chk_in_param(n_imgs,varargin{1});
        params.imag = chk_in_param(n_imgs,varargin{1});
    case {'dual_complex','dual_conformal'}
        params.type = 'dual_conformal';
        params.i{1} = chk_in_param(n_imgs,varargin{1})+ 1i.*chk_in_param(n_imgs,varargin{1});
        params.i{2} = chk_in_param(n_imgs,varargin{2})+ 1i.*chk_in_param(n_imgs,varargin{2});
    case 'Fourier'
        params.type = 'Fourier';
        params.fourier_coeffs = varargin{1};
        for k=1:1:size(params.fourier_coeffs,1)
            params.fourier_max_dist_X(k,:) = chk_in_param(n_imgs,varargin{2});
            params.fourier_max_dist_Y(k,:) = chk_in_param(n_imgs,varargin{end});
        end
end
end

function [out] = chk_in_param(n_imgs, in)
if max(size(in))==n_imgs
    out = in;
elseif all(size(in) == [1 2]) || all(size(in) == [2 1])
    out = rand(n_imgs,1) * (in(2)-in(1)) + in(1);
elseif numel(in)==1
    out = ones(n_imgs,1) .* in;
end
end

function [mdl] = conformal_mapping(mdl, fact_r, fact_i)
z_nodes = mdl.nodes(:,1) + 1i.*mdl.nodes(:,2);
z_mapped = z_nodes + (fact_r+1i*fact_i)*z_nodes.^2;
mdl.nodes(:,1) = real(z_mapped); mdl.nodes(:,2) = imag(z_mapped);
end

function [mdl] = dual_conformal_mapping(mdl, def1, def2)
z_nodes = mdl.nodes(:,1) + 1i.*mdl.nodes(:,2);
z_nodes = z_nodes + (def1)*z_nodes.^2;
z_nodes = z_nodes + (def2)*z_nodes.^2;
mdl.nodes(:,1) = real(z_nodes); mdl.nodes(:,2) = imag(z_nodes);
end

function [mdl] = distort_fourier(mdl,coefs)
geo.xy = coefs;
geo.z_mag = 1;
% Rotate -90 degrees
[theta, rho] = cart2pol(mdl.nodes(:,1), mdl.nodes(:,2));
[X, Y] = pol2cart(theta+deg2rad(-90), rho);
mdl.nodes(:,1) = X; mdl.nodes(:,2) = Y;
% Apply distortion
mdl = deform_cylinder(mdl,geo);
% Rotate 90 degrees
[theta, rho] = cart2pol(mdl.nodes(:,1), mdl.nodes(:,2));
[X, Y] = pol2cart(theta+deg2rad(90), rho);
mdl.nodes(:,1) = X; mdl.nodes(:,2) = Y;
end
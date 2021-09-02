function [ error_vec, mean_error ] = calc_error_MSE( target, reconstructed, varargin )
%CALC_ERRORMSE Compute the error, MSE method
%   Compute the error
%	Output the error as a vector, and the mean value of this vector
%	@Input :
%		target
%		reconstructed
%
%		varargin: additional arguments, can be one of the following strings
%			'area': compute the size of the elements first, useful if the
%				model has large and small elements
%			'volume': same as area, for 3D models
%			'normalize': In case the models are not normalized
%
%			'region': Equation specifying the error region
%				Next argument is the equation, string line or cell array
%				that will be used by EIDORS function elem_select
%			'list': List of elements on which we should calculate the error
%			'ellipse': creates an ellipse surrounding the target
%				Next parameter is a 5*1(2D) or 8*1(3D) array containing the parameters:
%					2 or 3 radii (cartesian coordinates)
%					1 or 2 Rotation angle of the ellipse (in degrees)
%					Position of the CoG of the ellipse (cartesian coordinates)
%			'half_ellipse': Half an ellipse, considers the part of the
%			ellipse close to the probe
%				Parameters are identical to 'ellipse'
%			'cylinder': Cylinder along the z-axis
%				Next parameter is a 5*1(2D) or 6*1(3D) array containing the parameters:
%					2 radii over the x and y axis (cartesian coordinates)
%					1 or 2 Rotation angle of the ellipse (in degrees)
%					Position of the CoG of the cylinder (x and y) (cartesian coordinates)
%			'aperture': Defines an aperture angle around the location of
%			the target
%				Next argument is the aperture angle, in degrees, aligned on
%				the x-axis
%			'function': Designed to use any drawing function to select the
%			ROI
%				Next argument is a cell array containing respectively:
%					- The drawing function (anonymous function)
%					- The CoG, size, and rotation angle, they can be
%					divided into 3 cells or sent as a 4*1 (2D) or 5*1 (3D) array
%			'weight': Followed by a percentage: Indicates the importance of
%				the specific region versus the elements out of this region
%				default: 1 (100%), elements out of this region are ignored
%
%	@Output :
%		error_vec: Vector of errors of each element of the mesh
%		mean_error: Actual MSE error, average of error vector error_vec
%

debug = false;
opts = process_args([],varargin{:});

% Basic calculation
error_vec = abs(get_img_data(reconstructed)-get_img_data(target)).^2; % ./ target.elem_data;
mean_error = mean(error_vec);

% Should be here to consider the other input arguments
if opts.normalize
    error_vec = error_vec / max(error_vec);
    mean_error = mean(error_vec);
end

if opts.calc_surf
    area_vec = abs(calc_elements_area(target.fwd_model));
    error_vec = (error_vec(:) .* (area_vec(:)))';
    mean_error = sum(error_vec)/(sum(area_vec));
else
    area_vec = ones(size( target.fwd_model.elems ));
end

selected = [];
if opts.use_region
    err_reg = opts.region_param.region;
    if ischar(err_reg); err_reg = inline(err_reg,'x','y','z'); end;
    selected = [selected, elem_select(target.fwd_model,err_reg)>0.5];
end
if opts.use_list
    tot_sz = size(target.fwd_model.elems,1);
    curr_sel = false(tot_sz,1);
    curr_sel( opts.list.sel_elems ) = true; %list of selected elements on which compute the error
    selected = [selected, curr_sel];
end
if opts.mk_ellipse
    n_dims = size(target.fwd_model.nodes,2);
    radius = opts.ellipse.params([1:n_dims]);
    if n_dims == 2
        ang_rot = repmat( opts.ellipse.params(3), 2,1);
        [ang_pos, center] = cart2pol(opts.ellipse.params(4), opts.ellipse.params(5));
        centerZ = center;
    elseif n_dims == 3
        ang_rot = repmat( opts.ellipse.params(4), 2,1);
        eidors_msg('3D Ellipse ROI: rotation angles are not all considered',3);
        [ang_pos, center] = cart2pol(opts.ellipse.params(6), opts.ellipse.params(7));
        centerZ=[center opts.ellipse.params(8)];
    end
    ang_pos = repmat( ang_pos, n_dims,1);
    [ell, ~] = createEllipseGeneral(1,target.fwd_model,1,radius,radius,center,centerZ,ang_rot,ang_pos,2);
    selected = [selected, (ell~=1)];
end
if opts.mk_function
    CoG = opts.use_fcn.CoG;
    sz = opts.use_fcn.sz;
    ang_rot = opts.use_fcn.ang_rot;
    [ang_pos, radius] = cart2pol(CoG(1),CoG(2));
    ang_pos = rad2deg(ang_pos);
    [ell, ~] = opts.use_fcn.fcn(1,target.fwd_model,1,sz,sz,radius,radius,ang_rot,ang_pos,2);
    selected = [selected, (ell~=1)];
end
if opts.mk_cylinder
    n_dims = size(target.fwd_model.nodes,2);
    radius = opts.cylinder.params([1:2]);
    ang_rot = repmat( opts.cylinder.params(3), 2,1);
    if n_dims == 2
        [ang_pos, center] = cart2pol(opts.cylinder.params(4), opts.cylinder.params(5));
    elseif n_dims == 3
        warning('Rotation angles are not all considered');
        [ang_pos, center] = cart2pol(opts.cylinder.params(5), opts.cylinder.params(6));
    end
    ang_pos = repmat( ang_pos, 2,1);
    [cyl, ~] = createCylinder(1,target.fwd_model,2,radius,radius,center,center,ang_rot,ang_pos,1);
    selected = [selected, (cyl~=2)];
end
if opts.mk_aperture
    ang = deg2rad( opts.aperture.angle/2 );
    err_reg{1} = inline(['y<x*',num2str(tan(ang))],'x','y','z');
    err_reg{2} = inline(['y>x*',num2str(tan(-ang))],'x','y','z');
    err_reg{3} = inline(['x>0'],'x','y','z');
    selected = [selected, elem_select(target.fwd_model,err_reg)>0.5];
end

if ~isempty(selected)
    selected = all(selected,2);
    err_in = sum(error_vec(selected))/(sum(area_vec(selected)));
    err_out = sum(error_vec(~selected))/(sum(area_vec(~selected)));
    mean_error = err_in * opts.weight + err_out * (1-opts.weight);
end

if debug
    dbgImg = mk_image(target.fwd_model,1);
    dbgImg.elem_data( selected ) = 2;
    figure; show_fem(dbgImg,1); title('Selected elements to determine error');
end

end

function args = process_args( args, varargin )
if isempty(args) %First call
    args.normalize = false;
    args.calc_surf = false;
    args.use_region = false;
    args.use_list = false;
    args.mk_ellipse = false;
    args.mk_cylinder = false;
    args.mk_aperture = false;
    args.weight = 1; % Don't change this to keep compatibility with old code
    args.mk_function = false;
end
n_params = length(varargin);
k=1;
while (k<=n_params)
    opt = varargin{k};
    if ischar(opt)
        switch opt
            case 'normalize'
                args.normalize = true; k=k+1;
            case 'area'
                args.calc_surf = true; k=k+1;
            case 'volume'
                args.calc_surf = true; k=k+1;
                
                % Different methods to select a parameter
            case 'region'
                args.use_region = true;
                args.region_param.region = varargin{k+1}; k=k+2;
            case 'list'
                args.use_list = true;
                args.list.sel_elems = varargin{k+1}; k=k+2;
            case 'ellipse'
                args.mk_ellipse = true;
                args.ellipse.params = varargin{k+1}; k=k+2;
            case 'half_ellipse'
                args.mk_ellipse = true;
                args.use_region = true;
                args.ellipse.params = varargin{k+1};
                args.region_param.region{1} = inline(['x<',num2str(varargin{k+1}(4))],'x','y','z');
                args.region_param.region{2} = inline('x>0','x','y','z');
                k=k+2;
            case 'cylinder'
                args.mk_cylinder = true;
                args.cylinder.params = varargin{k+1}; k=k+2;
            case 'aperture'
                warning('Aperture angle assumes the object is on the x-axis');
                args.mk_aperture = true;
                args.aperture.angle = varargin{k+1};
                if ((varargin{k+1}<0) || (varargin{k+1}>360))
                    error('Aperture angle should be a real number from 0 to 360');
                end
                k=k+2;
            case 'function'
                args.mk_function = true;
                if iscell(varargin{k+1})
                    args.use_fcn = get_params_target(varargin{k+1}); k=k+2;
                else
                    args.use_fcn = get_params_target(varargin{k+1},varargin{k+2}); k=k+3;
                end
                
                % Weight the target regions
            case 'weight'
                if ((varargin{k+1}<0) || (varargin{k+1}>1))
                    error('weight should be a real number from 0 to 1 (percentage)');
                end
                args.weight = varargin{k+1}; k=k+2;
        end
    elseif iscell(opt)
        args = process_args(args, opt{:}); % Here for backward compatibility of the code
        k=k+1;
    else
        eidors_msg(['Function ',mfilename,': Parameter not understood'],2);
        k=k+1;
    end
end
end

% If you use a specific drawing function to get the properties of the
% target, this function aims to obtain the correct parameters
function params =  get_params_target (varargin)
if (nargin==1) % Preferred method
    props = varargin{1};
    params.fcn = props{1};
    if iscell(props{2})
        params.CoG = props{2};
        params.sz = props{3};
        params.ang_rot = props{4};
        return;
    else
        lst = props{2};
    end
else
    params.fcn = varargin{1};
    lst = varargin{2};
end
if length(lst) == 4;
    params.CoG = lst(1:2);
    params.sz = lst(3);
    params.ang_rot = lst(4);
elseif length(lst)==5 % 3D
    params.CoG = lst(1:3);
    params.sz = lst(4);
    params.ang_rot = lst(5);
else
    error('Parameters not understood');
end
end

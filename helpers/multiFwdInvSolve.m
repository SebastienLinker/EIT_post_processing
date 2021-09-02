function [ invImgs, vi, vh ] = multiFwdInvSolve( model, homo_img, inhomo_imgs, num, bkgnd, parallel, voltages, do_inv )
%MULTIFWDINVSOLVE Perform forward and inverse solvers for multiple images
%   Then you don't have to write a loop in your main code
%	Input:
%		model: the EIDORS model
%		homo_img: Homogeneous image
%		inhomo_imgs: Inhomogenous images
%		num: number of inhomogeneous images
%	`	bkgnd: to adjust conductivity of the resulted images, set to 0 to
%		parallel: Use parallel computing (default false, broken on latest Matlab versions)
%		voltages: For real experiment, a N*2 array of N measurements: [vh; vi]
%       do_inv: Also solves inverse problem
%	Output:
%		invImgs: the result of the inverse solvers in a 1*num array
%		vi: an array containing the different measurements

if ~exist('parallel','var'); parallel = false; end
if ~exist('voltages','var'); voltages = []; end
if ~exist('do_inv','var'); do_inv = true; end

model.inv_solve.scale_solution.offset = bkgnd;
inv_in_loop = true;
% Pre calculate the reconstruction matrix
solver = model.solve;
if do_inv && ischar(solver)
    if strcmp(solver, 'eidors_default'); solver = eidors_default('get','inv_solve'); end
    solver = str2func(solver);
end
if do_inv
    calcRM = false;
    calcParamsPDIPM = false;
    if isequal(solver, @inv_solve_diff_GN_one_step)
        calcRM = true;
        model.solve = @inv_solve_diff_GN_one_step_no_RM;
        clean_mdl = true;
        readd = struct('solve',model.solve,'inv_solve',model.inv_solve);
        model = rmfield(model, {'solve','inv_solve'});
        inv_in_loop = false;
    elseif isequal(solver,@inv_solve_diff_pdipm)
        calcParamsPDIPM = true;
        model.solve = @inv_solve_diff_pdipm_fast;
        clean_mdl = false;
    elseif isequal(solver,@inv_solve_ANN)
        inv_in_loop = false;
    elseif isequal(solver,@inv_solve_post_proc)
        lin_solve = model.inv_solve.linear_solve;
        if ischar(lin_solve)
            if strcmp(lin_solve, 'eidors_default')
                lin_solve = eidors_default('get','inv_solve');
            end
            lin_solve = str2func(lin_solve);
        end
        if isequal(lin_solve,@inv_solve_diff_GN_one_step)
            calcRM = true;
            model.inv_solve.linear_solve = @inv_solve_diff_GN_one_step_no_RM;
            clean_mdl = true; % For eidors cache
            readd = struct('solve',model.solve,'inv_solve',model.inv_solve,...
                'inv_solve_post_proc',model.inv_solve_post_proc);
            model = rmfield(model, {'solve','inv_solve', 'inv_solve_post_proc'});
            inv_in_loop = false;
        elseif isequal(lin_solve,@inv_solve_diff_pdipm) || isequal(lin_solve,@inv_solve_diff_pdipm_fast)
            calcParamsPDIPM = true;
            model.inv_solve.linear_solve = @inv_solve_diff_pdipm_fast;
            clean_mdl = false;
        else
            clean_mdl = false;
        end
    end
    if calcRM
        eidors_msg('Pre-calculating the Reconstruction Matrix ...',2);
        RM = get_RM_GN_one_step( model );
        eidors_msg('RM calculated',2);
        if clean_mdl
            f_names = fieldnames(readd);
            for k=1:1:length(f_names)
                try model = setfield(model, f_names{k},readd.(f_names{k})); end
            end
        end
        model.inv_solve.RM = RM;
    elseif calcParamsPDIPM
        eidors_msg('Pre-calculating the parameters for PDIPM method ...',2);
        [model.inv_solve.J, model.inv_solve.L, model.inv_solve.LTL, model.inv_solve.JTW] = ...
            get_J_L_LTL_JTW( model );
        eidors_msg('Parameters calculated',2);
    end
end

% bound_homo = fwd_solve( homo_img );
do_fwd = isempty(voltages);
do_homo = true;
try
    if do_inv && ~do_fwd; [voltages(:,1).meas]; end % Outputs error
catch; inv_in_loop = true;
end
if isfield(model,'do_homo') && model.do_homo==false; do_homo=false; end
if isfield(model,'reconst_type') && strcmp(model.reconst_type,'static'); do_homo=false; end
if ~do_fwd
    if do_homo
        [vi, vh] = reshape_voltages(voltages, do_homo);
    else
        [vi] = reshape_voltages(voltages, do_homo);
    end
else
    vh = repmat( struct('meas',[],'time',NaN,'name','','type','data'), 1,num);
    vi = repmat( struct('meas',[],'time',NaN,'name','','type','data'), 1,num);
end
if do_inv
    invImgs = repmat( struct('current_params',[],'elem_data',[],'fwd_model',[],'info',[],...
        'name','','type',''), 1,num);
end

log_lvl = eidors_msg('log_level');
eidors_msg('log_level',1); % Then we don't print out tons of stuffs

if (license('test','distrib_computing_toolbox') && parallel && (num>50))
    parallelpool('open');
    parfor k = 1:1:num
        % Forward
        if do_fwd
            [vh_tmp, vi_tmp] = calc_loop_fwd (homo_img, inhomo_imgs(k), do_homo);
            vi(k) = vi_tmp;
            if do_homo; vh(k) = vh_tmp; end;
        else
            vh_tmp = vh(k); vi_tmp = vi(k);
        end
        % Inverse
        if do_inv && inv_in_loop
            sol_tmp = calc_loop_inv(model, vh_tmp, vi_tmp);
            invImgs(k) = sol_tmp;
        end
    end
    %     if (matlabpool('size')); matlabpool close; end
else
    for k = 1:1:num
        % Forward
        if do_fwd
            [vh_tmp, vi_tmp] = calc_loop_fwd (homo_img, inhomo_imgs(k), do_homo);
            if do_homo; vh(k) = vh_tmp; end;
            vi(k) = vi_tmp;
        else
            vi_tmp = vi(k);
            if do_homo; vh_tmp = vh(k); else vh_tmp=[]; end;
        end
        % Inverse
        if do_inv && inv_in_loop
            sol_tmp = calc_loop_inv(model, vh_tmp, vi_tmp);
            invImgs(k) = sol_tmp;
        end
    end
end

% Solve inverse problems outside the loop
% Much faster but doesn't work well with all solvers
if do_inv && ~inv_in_loop
    if do_homo;	sol_tmp = calc_loop_inv(model, vh, vi);
    else sol_tmp = calc_loop_inv(model, vh); end
    invImgs = repmat( sol_tmp, 1,length(vh) );
    for k=1:1:length(vh); invImgs(k).elem_data = sol_tmp.elem_data(:,k); end
end

if ~do_inv; invImgs = []; end;
if ~exist('vh','var'); vh = []; end;

eidors_msg('log_level',log_lvl);
if isfield(homo_img,'calc_colours')
    for k = 1:1:num
        invImgs(k).calc_colours = homo_img.calc_colours;
    end
end

end

function [vh, vi] = calc_loop_fwd (homo, inh, do_homo)
if do_homo
    vh_t = fwd_solve( homo ); % Better to do this here, in case we add random noise
end
vi_t = fwd_solve( inh );
% Just because of EIDORS 3.8 and additional fields
vh.meas=vh_t.meas; vh.time=vh_t.time; vh.name=vh_t.name; vh.type=vh_t.type;
vi.meas=vi_t.meas; vi.time=vi_t.time; vi.name=vi_t.name; vi.type=vi_t.type;
end

function [inv_sol] = calc_loop_inv(mdl, vh, vi)
% Inverse problem
if strcmp(mdl.reconst_type, 'difference');
    tmp = inv_solve( mdl, vh, vi );
else
    tmp = inv_solve( mdl, vi );
end

inv_sol.current_params = tmp.current_params;
inv_sol.elem_data = tmp.elem_data;
inv_sol.fwd_model = tmp.fwd_model;
inv_sol.info = tmp.info;
inv_sol.name = tmp.name;
inv_sol.type = tmp.type;
end

% Reshape vh and vi arrays, in case they are provided
function [vi, vh] = reshape_voltages(voltages, do_homo)
if isstruct(voltages)
    for k = 1:1:size(voltages,1)
        v = voltages(k,:);
        if do_homo
            vh(k) = v(1); vi(k) = v(2);
        else vi(k) = v;
        end
    end
else
    v=voltages;
    % 		warning('Please verify this case');
    if do_homo
        vh.meas = v(:,1); vh.type = 'data';
        vi.meas = v(:,2); vi.type = 'data';
    else
        vi.meas = v(:,1); vi.type = 'data';
    end
end
end

function pool = parallelpool( action, pool )
v = ver('MATLAB');
v = str2double(v.Version);
if (v>=8.2)
    use_parpool = true;
else
    use_parpool = false;
end
switch action
    case 'open'
        if parallelpool('test'); return; end
        if use_parpool; pool = parpool('local',4);
        else; matlabpool local 4; pool = [];
        end
    case 'close'
        if ~parallelpool('test'); return; end
        if use_parpool; delete(gcp);
        else; matlabpool close; pool = [];
        end
    case 'test'
        if use_parpool;
            tmp = gcp('nocreate');
            if ~isempty(tmp); pool = true; else pool = false; end
        else
            if ~matlabpool('size'); pool = true; else pool = false; end
        end
    otherwise
        error('Unknown action');
end
end
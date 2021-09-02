function [ inh, inh_n, solve, solve_n, vh, vi ] = solveAndNormalize( conduct, mdl, bk, parallel, voltages )
%SOLVEANDNORMALIZE Helper function which intents to create the models, solve
%and normalize the images
%   Purpose: Avoid multiple copy-pastes
%
%	@input conduct: conductivity of each element for each image
%	@input imdl: inverse model (includes forward)
%		note: If it is a forward model, then only the forward problems are
%		solved
%	@input bkgnd: homogeneous conductivity (background)
%	@input parallel: solve the inverse problems in parallel
%	@input solver: inverse solver to use
%	@input prior: prior probability to use
%	@input voltages: For real experiment, a N*2 array of N measurements: [vh; vi]
%
%	@output inhomo, inhomo_norm: original image, and normalization
%	@output solve, solve_norm: inverse solution, and normalization
%	@output dv: voltage difference
%

if ~exist('voltages','var'); voltages = []; end

if isfield(mdl,'fwd_model')
    fmdl = mdl.fwd_model; do_inv = true;
    imdl = mdl;
else
    fmdl = mdl; do_inv = false;
end
img_homo = mk_image( fmdl, bk);
img_homo.current_params = []; img_homo.info = [];

if numel(voltages)
    if isstruct(voltages)
        n_images = size(voltages,1);
    elseif ismatrix(voltages)
        n_images = 1;
    else
        error('Input voltages not understood');
    end
else
    n_images = numel (conduct(1,:));
end

for k=1:1:n_images
    inh(k) = img_homo;
    if isscalar(conduct)
        tmp = mk_image(img_homo.fwd_model, conduct); % In case conducts is a scalar
        tmp.current_params = []; tmp.info = [];
        inh(k) = tmp;
    else
        inh(k).elem_data = conduct(:,k);
    end
    inh(k).fwd_model.normalize_measurements = 1;
    if isempty(voltages)
        inh_n(k) = normalize_mdl( inh(k) );
    else; inh_n = []; % Just to make sure it won't crash at this point
    end
end

if do_inv; imdl.inv_solve.scale_solution.offset = bk; end
% Check system matrix: Using EIDORS 3.7.1: Inverse solver assumes 1st
% order matrix, but can have higher order forward solvers
fwd_sol = fmdl.solve;
if isa(fwd_sol,'function_handle')
    fwd_sol = func2str(fwd_sol);
end
if strcmp(fwd_sol, 'eidors_default');
    fwd_sol = eidors_default('get','fwd_solve');
end
if do_inv && strcmp(fwd_sol, 'fwd_solve_higher_order')
    eidors_msg(['You selected higher order forward solver, but EIDORS v3.7 ', ...
        'assumes a 1st order matrix in inverse solver. To avoid errors, 1st ', ...
        'order forward solver is used'],1);
    mdl.fwd_model.solve = @fwd_solve_1st_order;
    mdl.fwd_model.system_mat = @system_mat_1st_order;
    mdl.fwd_model.jacobian = @jacobian_adjoint;
end

if do_inv; mdl.do_homo = ~(strcmp(mdl.reconst_type,'absolute')); end
[solve, vi, vh] = multiFwdInvSolve(mdl, img_homo, inh, n_images, bk, parallel, voltages, do_inv);
if do_inv
    try
        solve_n = normalize_mdl( solve, bk );
    catch
        solve_n=[];
    end
else
    solve_n = [];
end

end
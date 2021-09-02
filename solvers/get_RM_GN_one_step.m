function [ RM ] = get_RM_GN_one_step( imdl )
%GET_RM_GN_ONE_STEP Get reconstruction matrix for difference inverse
%solver inv_solve_diff_GN_one_step_no_RM
%   Input: imdl, the inverse model
%	Output: RM: The reconstruciton matrix

imdl = prepare_model( imdl );
RM = eidors_cache(@get_RM, imdl,'get_RM_GN_one_step' );

end

% Obtained from inv_solve file
function mdl = prepare_model( mdl )
fmdl = mdl.fwd_model;
fmdl = mdl_normalize(fmdl,mdl_normalize(fmdl));
if ~isfield(fmdl,'elems')
    return;
end

fmdl.elems  = double(fmdl.elems);
fmdl.n_elem = size(fmdl.elems,1);
fmdl.n_node = size(fmdl.nodes,1);
if isfield(fmdl,'electrode')
    fmdl.n_elec = length(fmdl.electrode);
else
    fmdl.n_elec = 0;
end

% So we can use cached values, noise doesn't matter here
if (isa(fmdl.solve,'function_handle') && strcmp(func2str(fmdl.solve),'fwd_solve_add_noise')) || ...
        (ischar(fmdl.solve) && strcmp(fmdl.solve,'fwd_solve_add_noise'))
    try fmdl.solve = fmdl.add_noise.solve;
    catch; fmdl.solve = 'eidors_default'; end
    if isfield(fmdl,'add_noise'); fmdl = rmfield(fmdl,'add_noise'); end
end

mdl.fwd_model= fmdl;
if ~isfield(mdl,'reconst_type')
    mdl.reconst_type= 'difference';
end

end

% Inspired from file inv_solve_diff_GN_one_step
% No need to check the cache, it is assumed to be too large anyway
function RM = get_RM( inv_model )
% The one_step reconstruction matrix is cached
%    RM = eidors_obj('get-cache', inv_model, 'inv_solve_diff_GN_one_step');
%    if ~isempty(RM)
%        eidors_msg('inv_solve_diff_GN_one_step: using cached value', 3);
%        return;
%    end

img_bkgnd= calc_jacobian_bkgnd( inv_model );
img_bkgnd.fwd_model.jacobian_SbS.jacobian_func = img_bkgnd.fwd_model.jacobian;
img_bkgnd.fwd_model.jacobian = @jacobian_step_by_step;
J = calc_jacobian( img_bkgnd );
J = single(J);

RtR = calc_RtR_prior( inv_model );
W   = calc_meas_icov( inv_model );
hp  = calc_hyperparameter( inv_model );

% This line needs too much memory, slow down the calculation
%    RM= (J'*W*J +  hp^2*RtR)\J'*W;
% This code does the same, needs more operations but needs less memory, then it is faster
JT = J';
if isequal(W,speye(size(W,1)))
    % Output of default inverse covariance matrix, so it is equivalent
    % to multiplying it by identity matrix (i.e. does nothing)
    tmp = JT*J;
else
    tmp = JT*full(W); tmp=tmp*J; % This one might be slow
end
clear J;
hp2RtR = hp^2*RtR; clear hp RtR;
if issparse(hp2RtR) && ~isa(tmp,'single') % prior_exponential_covar_3D outputs nonsparse matrix
    % Take advantage of sparsity, don't convert to full matrix
    ind = (hp2RtR~=0); % todo: try nonzeros(), maybe it is faster
    tmp(ind) = tmp(ind) + hp2RtR(ind);
    clear ind;
elseif isa(tmp,'single')
    ind = (hp2RtR~=0);
    tmp(ind) = tmp(ind) + single(full(hp2RtR(ind)));
    clear ind;
else
    tmp = tmp + hp2RtR;
end
clear hp2RtR; % sparse matrix, but we don't need it anymore

% 	curr_RM = tmp\JT; clear tmp JT; % Still the slowest line (need splitting)
max_meas = 50;
n_loops = ceil(size(JT,2)/max_meas);
JTorRM = JT; clear JT; % Keep the same space in memory
for k = 1:1:n_loops
    curr_max_pos = min(k*max_meas, size(JTorRM,2));
    curr_pos = (k-1)*max_meas+1:1:curr_max_pos;
    JTorRM(:,curr_pos) = tmp\JTorRM(:,curr_pos);
end
curr_RM = JTorRM;
clear JTorRM tmp;

if ~isequal(W,speye(size(W,1)))
    curr_RM = curr_RM*full(W);
end
RM = curr_RM; clear curr_RM;

%    eidors_obj('set-cache', inv_model, 'inv_solve_diff_GN_one_step', RM);
%    eidors_msg('inv_solve_diff_GN_one_step: setting cached value', 3);
end
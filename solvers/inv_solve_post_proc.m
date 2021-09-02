function img= inv_solve_post_proc( inv_model, varargin)
% INV_SOLVE_DIFF_POST_PROC Post processing with ANN
% img= inv_solve_diff_post_proc( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%
% Please set these parameters in the inverse model
% inv_model.inv_solve.linear_solve: Linear inverse solver
% inv_model.inv_solve.ann: the neural network for post processing
% inv_model.inv_solve.normalize: boolean, should we normalize before PP? (default false)

% (C) 2015 Sebastien Martin

solver = inv_model.inv_solve.linear_solve;
pp_solver = inv_model.solve;

bk = inv_model.inv_solve.scale_solution.offset;
inv_model.solve = solver;

inv_solve_post_proc = inv_model.inv_solve_post_proc;

inv_model_lin = rmfield(inv_model,'inv_solve_post_proc');
if isfield(inv_solve_post_proc,'dual_mesh')
    inv_model_lin.fwd_model = inv_solve_post_proc.dual_mesh;
    inv_model_lin.fwd_model.stimulation = inv_model.fwd_model.stimulation;
    try inv_model_lin.fwd_model.meas_select = inv_model.fwd_model.meas_select; end
    inv_model_lin.inv_solve.select_parameters = [1:1:size(inv_model_lin.fwd_model.elems,1)];
    try inv_model_lin.inv_solve = rmfield( inv_model_lin.inv_solve,'RM' ); end
    try inv_model_lin.inv_solve = rmfield( inv_model_lin.inv_solve,'ann_post_proc' ); end
    try inv_model_lin.inv_solve = rmfield( inv_model_lin.inv_solve,'linear_solve' ); end
    warning('off','EIDORS:lateRM');
end
img = inv_solve(inv_model_lin, varargin{:});
if isfield(inv_solve_post_proc,'dual_mesh'); warning('on','EIDORS:lateRM'); end

inv_model.solve = pp_solver;

try ann = inv_model.inv_solve.ann_post_proc; catch; ann = inv_model.parameters.ann_post_proc; end

% Should we normalize the image before aplying ANN?
try
    if inv_model.inv_solve.normalize;
        img = normalize_mdl(img,bk);
    end;
catch
end

if isfield(inv_solve_post_proc,'dual_mesh')
    img.fwd_model = inv_model.fwd_model;
end

% It should not be at the expected position fmdl.coarse2fine, since we want
% to use the large FEM for linear solver
if isfield(inv_solve_post_proc,'c2f')
    c2f = inv_solve_post_proc.c2f;
    if size(c2f,2)==size(img.elem_data,1)
        img.elem_data = c2f*img.elem_data;
    else eidors_msg('@@@ Coarse to fine mapping: Size of matrices are not consistent, mapping ignored');
    end
end

img.name= [img.name,' and post-processing'];
if isfield(inv_model.inv_solve,'switch_range') && inv_model.inv_solve.switch_range
    img.elem_data = 1-img.elem_data;
end

try do_nodes = inv_solve_post_proc.do_nodes;
catch; do_nodes = false; end; %Backward compatibility
if do_nodes
    img.node_data = ann(inv_solve_post_proc.e2n * img.elem_data);
    img.elem_data = inv_solve_post_proc.n2e * img.node_data;
else
    img.elem_data = ann(img.elem_data);
end
if isfield(inv_model.inv_solve,'switch_range') && inv_model.inv_solve.switch_range
    img.elem_data = 1-img.elem_data;
end

end
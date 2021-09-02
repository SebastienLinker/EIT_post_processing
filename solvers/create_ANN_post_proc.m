function [ imdl, net ] = create_ANN_post_proc( imdl, img_solve, img_inhomo, net_params, bk, varargin )
%CREATE_ANN_POST_PROC Create ANN for postprocessing
%	Returns the inverse model ready for post processing and the ANN itself
%   Requires MATLAB's neural network toolbox
%
% (C) 2015/05/23 Sebastien Martin

opts = parse_varargin(varargin{:});
net_params = configure_ANN(net_params);

n_elements = length(img_inhomo(1).elem_data);
n_images = length(img_inhomo);

if opts.normalize
    img_solve = normalize_mdl(img_solve, bk);
    img_inhomo = normalize_mdl(img_inhomo, bk);
end

T = [img_inhomo(:).elem_data]; P = [img_solve(:).elem_data];

if opts.do_c2f
    try
        c2f= img_inhomo(1).fwd_model.coarse2fine;
        T = c2f * T;
    catch; eidors_msg('@@@ Coarse to fine mapping failed', 3);
    end
    try
        c2f= img_solve(1).fwd_model.coarse2fine;
        P = c2f * P;
    catch; eidors_msg('@@@ Coarse to fine mapping failed', 3);
    end
end

if opts.solve_nodes
    if opts.do_c2f
        [e2n_lin, n2e_lin] = node_elem_mapper(img_solve(1).fwd_model);
        [e2n_sol, n2e_sol] = node_elem_mapper(img_inhomo(1).fwd_model);
    else
        if size(P,1)==size(imdl.fwd_model.elems,1)
            [e2n_lin, n2e_lin] = node_elem_mapper(imdl.fwd_model);
        else
            [e2n_lin, n2e_lin] = node_elem_mapper(img_solve(1).fwd_model);
        end
        [e2n_sol, n2e_sol] = node_elem_mapper(img_inhomo(1).fwd_model);
    end
    P = e2n_lin*P;
    %backward compatibility? not sure we can use what's in 'catch' only
    try
        T = e2n_lin * T;
        use_n2esol = false;
    catch
        T = e2n_sol * T;
        use_n2esol = true;
    end
end

% Create neural networks
net = fitnet(net_params.hid_lay);
net.trainFcn = net_params.train;
net.layers{1}.transferFcn = net_params.trans{1};
net.layers{2}.transferFcn = net_params.trans{2};
net.trainParam.goal = net_params.goal;
net.trainParam.epochs = net_params.epochs;
net.trainParam.max_fail = net_params.max_fail;
net.trainParam.min_grad = net_params.min_grad;

[net,~] = train(net, P, T);

lin_solve = imdl.solve;
if ischar(lin_solve)
    if strcmp(lin_solve, 'EIDORS_default')
        lin_solve = eidors_default('get','inv_solve');
    end
    lin_solve = str2func(lin_solve);
end
if opts.getRM && isequal(lin_solve, @inv_solve_diff_GN_one_step)
    imdl.solve = @inv_solve_diff_GN_one_step_no_RM;
    eidors_msg('Pre-calculating the reconstruction matrix ...',2);
    imdl.inv_solve.RM = get_RM_GN_one_step(imdl);
    eidors_msg('Done',2);
end

imdl.inv_solve.ann_post_proc = net;
imdl.inv_solve.linear_solve = imdl.solve;
imdl.inv_solve.normalize = opts.normalize;
imdl.solve = @inv_solve_post_proc;
imdl.inv_solve_post_proc.do_nodes = opts.solve_nodes;
if opts.solve_nodes
    imdl.inv_solve_post_proc.e2n = e2n_lin;
    if ~use_n2esol
        imdl.inv_solve_post_proc.n2e = n2e_lin;
    else
        imdl.inv_solve_post_proc.n2e = n2e_sol;
    end
end

if opts.do_c2f
    imdl.inv_solve_post_proc.c2f = c2f;
else
    try
        imdl.inv_solve_post_proc = rmfield(imdl.inv_solve_post_proc,'c2f');
    catch
    end
end

if opts.dual_mesh.do
    imdl.inv_solve_post_proc.dual_mesh = opts.dual_mesh.mesh;
end

end

%		net_p: parameters of the ANN
%			-hid_lay: number of neurons in the hidden layer (default 10)
%			-train: training function (default traingdx)
%			-trans: A 2*1 string cell array with the name of the transfer
%			function (default {'tansig','tansig'})
%			-goal: (default 0.000001)
%			-epochs: (default 10000)
%			-max_fail: (default 500)
function net_p = configure_ANN(net_p)
if ~isfield(net_p,'hid_lay'); net_p.hid_lay=100; end
if ~isfield(net_p,'train'); net_p.train='traingdx'; end
if ~isfield(net_p,'trans'); net_p.trans={'tansig','tansig'}; end
if ~isfield(net_p,'goal'); net_p.goal=0.000001; end
if ~isfield(net_p,'epochs'); net_p.epochs=10000; end
if ~isfield(net_p,'max_fail'); net_p.max_fail=500; end
if ~isfield(net_p,'min_grad'); net_p.min_grad=1e-5; end
if ~isfield(net_p,'difference'); net_p.difference=false; end
end

function opts = parse_varargin(varargin)
%default parameters
opts.getRM = false;
opts.normalize = true;
opts.solve_nodes = false;
opts.do_c2f = false;
opts.dual_mesh.do = false;

k=1;
while k <= length(varargin)
    opt = varargin{k};
    switch opt
        case 'getRM'
            opts.getRM = true;
        case 'noRM'
            opts.getRM = false;
        case 'norm'
            opts.normalize = true;
        case 'no_norm'
            opts.normalize = false;
        case 'do_nodes'
            opts.solve_nodes = true;
        case 'do_c2f'
            opts.do_c2f = true;
        case 'dual_mesh'
            opts.dual_mesh.do = true;
            opts.dual_mesh.mesh = varargin{k+1}; k=k+1;
    end
    k = k+1;
end
end
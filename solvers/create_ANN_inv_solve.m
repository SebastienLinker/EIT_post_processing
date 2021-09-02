function [ imdl, net ] = create_ANN_inv_solve( imdl, voltages, img_inhomo, net_params, varargin )
%CREATE_ANN_INV_SOLVE Creates ANN as inverse solver
%	Returns the inverse model ready for solving inverse problem, and the ANN itself
%   Requires MATLAB's neural network toolbox
%
%	(C) 2015/06/22 Sebastien Martin

opts = parse_varargin(varargin{:});
net_params = configure_ANN(net_params);

n_elements = length(img_inhomo(1).elem_data);
n_images = length(img_inhomo);

T = [img_inhomo(:).elem_data];
try
    c2f= img_inhomo(1).fwd_model.coarse2fine;
    T = c2f * T;
catch
end

if isstruct(voltages)
    if size(voltages,2)==2;
        net_params.difference = true;
        P = calc_difference_data(voltages(:,1), voltages(:,2), imdl.fwd_model);
    else
        P = [voltages(:,1).meas];
    end
else
    P = voltages;
end

if opts.solve_nodes
    [e2n, n2e] = node_elem_mapper(img_inhomo(1).fwd_model);
    if size(e2n,2)~=size(T,1); [e2n, n2e] = node_elem_mapper(imdl.fwd_model); end
    T = e2n*T;
end

% Create neural networks
net = feedforwardnet(net_params.hid_lay, net_params.train);
net.trainFcn = net_params.train;
net.layers{1}.transferFcn = net_params.trans{1};
net.layers{2}.transferFcn = net_params.trans{2};
net.trainParam.goal = net_params.goal;
net.trainParam.epochs = net_params.epochs;
net.trainParam.max_fail = net_params.max_fail;
net.trainParam.min_grad = net_params.min_grad;

[net,~] = train(net, P, T);

imdl.inv_solve.ann = net;
imdl.reconst_type = 'difference';
if ~net_params.difference; imdl.reconst_type='absolute'; end
imdl.inv_solve.normalize = opts.normalize;
imdl.solve = @inv_solve_ANN;
imdl.inv_solve_ANN.do_nodes = opts.solve_nodes;
if opts.solve_nodes
    imdl.inv_solve_ANN.n2e = n2e;
end

end

%		net_p: parameters of the ANN
%			-hid_lay: number of neurons in the hidden layer (default 100)
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

k=1;
while k <= length(varargin)
    opt = varargin{k};
    switch opt
        case 'do_nodes'
            opts.solve_nodes = true;
    end
    k = k+1;
end
end

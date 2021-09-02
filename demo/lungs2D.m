%% Post processing 2D EIT bounded problems
% This file intents to train and test ANN for post processing real lungs
% data
% Here we do not make use of function newrb or newrbe, known for their poor
% generalization capabilities
%
% (C) 2015/10/14 Sebastien Martin

clearvars; % close all;

%% Parameters
do_simulation = true;
% Tank
n_elecs = 16; % # of electrodes
elec_diam = 0.4; % diameter of one electrode (mm)
bk = 0.0001;

% Targets
n_imgs = 1000;
n_imgs_test = 5;
max_sz= [0.6];
min_sz= [0.2];
min_pos = 0.1; max_pos = 0.7;% Distance from center
cond_r = [0.8 1]; rot_angle = [0 180]; pos_angle = {[-90 90],[90 270]};
fcn_mk_targ = @createEllipseGeneral; % For training ANN
n_targs = 2; fcn_targ_test = @createEllipseGeneral; overlap = false; % Testing

% Stimulation
opts_stim = {'no_meas_current','no_rotate_meas','do_redundant'};
amp = -0.1; %mA
noisy = false;

% Mesh distortion
distort.large_dist = false; % Large distortions
distort.simple = true; % Generate one simple model of the lungs
distort.train = true;
dist_type = 'Fourier'; %'complex'; % 'ellipse';
dist_params.X = [0.9 1.1]; dist_params.Y = [1 1.05];
dist_params.Fourier_C = shape_library('get','adult_male','boundary');

% Inverse solver
prior = @prior_movement;
jacobian = @jacobian_movement;
do_parallel = false;

% ANN
createNN = true; % True if you want to use a NN, complete
ann_param.fcn = 'traingdx';
ann_param.goal = 0.0001;
ann_param.min_grad = 1e-10;
ann_param.epochs = 10000;
ann_param.trans = {'radbas','tansig'};
ann_param.max_fail = 500; % # of consecutive iterations without improved performance
ann_param.hid_lay = 400; %hidden layer
ann_param_nopp = ann_param;
ann_param_nopp.trans = {'radbas','tansig'};
ann_param_nopp.hid_lay = 400;

%% 2 Create tank model and inverse model
maxh = 0.05;
imdl_fwd = mk_common_model('f2C',n_elecs); % Keep a copy of the original one, in case of distortion
imdl_fwd.fwd_model = mdl_normalize(imdl_fwd.fwd_model,1);
fmdl_fwd = imdl_fwd.fwd_model;
fmdl_fwd.jacobian_bkgnd.value = bk;
fmdl_fwd.jacobian = jacobian;

imdl = mk_common_model('e2C',n_elecs);
imdl.solve = @inv_solve_diff_GN_one_step; %imdl.hyperparameter.value = 100;
imdl.fwd_model = mdl_normalize(imdl.fwd_model,1);
imdl.RtR_prior = prior;
fmdl = imdl.fwd_model;
fmdl.jacobian_bkgnd.value = bk;
fmdl.jacobian = jacobian;
fmdl_bak = fmdl; %Backup the model
c2f_map = mk_coarse_fine_mapping(fmdl,fmdl_fwd);
fmdl_fwd.coarse2fine = c2f_map;

n_elements = size(fmdl_fwd.elems,1);
% Forward models
if noisy
    fmdl_fwd.solve = @fwd_solve_add_noise;
    fmdl_fwd.add_noise.filter = true;
    fmdl_fwd.parameters.isLungs = true;
    % 	fmdl_fwd.add_noise.noise = get_draeger_noise;
    % 	fmdl.add_noise.solve = @fwd_solve_higher_order;
    % 	fmdl.approx_type = 'tri6';
    % 	fmdl.system_mat= @system_mat_higher_order;
    % 	fmdl.jacobian = @jacobian_adjoint_higher_order;
end
% Inverse models
imdl.fwd_model = fmdl;
imdl.prior_movement.RegC.func = @prior_gaussian_HPF;
imdl.prior_movement.parameters = 0.0032;
imdl.inv_solve.select_parameters = [1:1:size(imdl.fwd_model.elems,1)];

if true
    switch dist_type
        case 'ellipse'
            args_dist = {dist_params.X, dist_params.Y};
        case 'complex'
            args_dist = {dist_params.conformal};
        case 'dual_complex'
            args_dist = dist_params.dual_conf(:);
        case 'Fourier'
            args_dist = {dist_params.Fourier_C, dist_params.X, dist_params.Y};
    end
end

%% Create ANNs
if createNN
    %% 3 Randomly creates images
    [ conds_NN, area_NN ] = createMultiTarg( n_imgs, n_targs, 0, fmdl_fwd, bk, @createEllipseGeneral,...
        max_sz,min_sz, min_pos,max_pos, rot_angle,pos_angle, cond_r);
    % 	[ conds_NN, area_NN ] = createLungs_train2D( n_imgs, fmdl_fwd, bk, cond_r );
    
    %% 3.5 Get the corresponding training data
    [stim,msel]= mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
    fmdl_fwd.stimulation = stim; fmdl_fwd.meas_select = msel;
    % Forward and inverse solvers
    if distort.train
        fmdl_fwd = mk_distortion(fmdl_fwd, n_imgs, dist_type, args_dist{:});
        for l = 1:1:n_imgs
            [inh_NN(l),~,~,~,t_vh(1,l),t_vi(1,l)] = solveAndNormalize(conds_NN(:,l), fmdl_fwd(l), bk, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    elseif (distort.large_dist || distort.simple)
        fmdl_fwd = mk_distortion(fmdl_fwd,1, dist_type, args_dist{1}, 1);
        [inh_NN,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_NN, fmdl_fwd, bk, do_parallel);
    else % Circular mesh, no distortion
        [inh_NN,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_NN, fmdl_fwd, bk, do_parallel);
    end
    if (distort.large_dist || distort.simple)
        imdl.fwd_model = mk_distortion(fmdl,1, dist_type, args_dist{1}, 1);
    end
    imdl.fwd_model.stimulation = stim; imdl.fwd_model.meas_select = msel;
    [~,~,t_sol,~,~,~] = solveAndNormalize(0,imdl,bk,do_parallel,[t_vh;t_vi]');
    
    %% 4 Actually creates 1 ANN for PP + 1 ANN with the classic method
    imdl_pp = create_ANN_post_proc(imdl, t_sol, inh_NN, ann_param, bk, 'no_norm', 'do_nodes', 'do_c2f');
    imdl_ann = create_ANN_inv_solve(imdl, [t_vh;t_vi]', inh_NN, ann_param_nopp, 'do_nodes');
    % Keep voltages for later use?
    vh_NN=t_vh; vi_NN=t_vi; sol_NN=t_sol; clear t_vh t_vi t_sol stim msel;
    disp(['Trained ANNs']);
    
    %% Clean up structures (and save lots of memory)
    f_to_rem = {'stimulation','meas_select'}; %fields to remove
    fmdl = rmfield(fmdl,f_to_rem); fmdl_fwd = rmfield(fmdl_fwd,f_to_rem);
    imdl.fwd_model = rmfield(imdl.fwd_model, f_to_rem);
    imdl_pp.fwd_model = rmfield(imdl_pp.fwd_model, f_to_rem);
    imdl_ann.fwd_model = rmfield(imdl_ann.fwd_model, f_to_rem);
    [inh_NN(:).fwd_model] = deal([]);
else
    % LOAD YOUR PRE TRAINED MODELS HERE
end

%% 5 Randomly generate new images
fmdl_fwd = fmdl_fwd(1); fmdl_fwd.nodes = imdl_fwd.fwd_model.nodes;
[conds_test, area_test] = createLungs2D( n_imgs_test, fmdl_fwd, bk, cond_r);
xyzr = [area_test(:,4) area_test(:,5) zeros(n_imgs_test,1) area_test(:,1)];
%% Shape deformation
if distort.simple
    fmdl_fwd = mk_distortion(fmdl_fwd, 1, dist_type, args_dist{1}, 1);
    fmdl = mk_distortion(fmdl_bak, 1, dist_type, args_dist{1}, 1);
elseif distort.large_dist
    fmdl_fwd = mk_distortion(fmdl_fwd, n_imgs_test, dist_type, args_dist{:});
    fmdl = mk_distortion(fmdl_bak, n_imgs_test, dist_type, args_dist{:});
else
    fmdl_fwd = repmat(fmdl_fwd,n_imgs_test,1);
    fmdl = repmat(fmdl_bak,n_imgs_test,1);
end

%% 5.5 If phantom measurements, load data
if ~do_simulation
    % LOAD YOUR MEASUREMENT DATA HERE
end

%% Make nonlinear model
imdl_nonl = imdl;
imdl_nonl.solve = @inv_solve_diff_pdipm;
imdl_nonl.parameters.term_tolerance= 0;
imdl_nonl.parameters.max_iterations = 5;
imdl_nonl.inv_solve_diff_pdipm.norm_data = 1;

%% 6 Apply ANN for post processing and linear solver (comparison purpose)
clear vh_test vi_test img_lin img_nonl img_ann img_pp;
[stim, msel] = mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
fmdl_fwd.stimulation = stim;
fmdl_fwd.meas_select = msel;
fmdl.stimulation = stim;
fmdl.meas_select = msel;
imdl.fwd_model=fmdl; imdl_nonl.fwd_model=fmdl;
imdl_pp.fwd_model=fmdl;
imdl_ann.fwd_model = fmdl;
% Get voltages first, in case we add a random noise, make sure the voltages are always the same
if do_simulation
    if distort.large_dist
        for l = 1:1:n_imgs_test
            [inh_t(l),~,~,~,t_vh(l),t_vi(l)] = solveAndNormalize(conds_test(:,l), fmdl_fwd(l), bk, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    elseif distort.simple
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl_fwd, bk, do_parallel);
    else
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl_fwd, bk, do_parallel);
    end
else
    t_vh = vh_ph; t_vi = vi_ph;
end
% Linear solver, Post-processing, classic ANN
[~,~,t_lin,~,~,~] = solveAndNormalize(0,imdl,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_ann,~,~,~] = solveAndNormalize(0,imdl_ann,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_pp,~,~,~] = solveAndNormalize(0,imdl_pp,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_nonl,~,~,~] = solveAndNormalize(0,imdl_nonl,bk,do_parallel,[t_vh;t_vi]');
% Free space
img_lin=t_lin;
img_pp=t_pp;
img_ann=t_ann;
img_nonl=t_nonl;
clear t_vh t_vi t_lin t_pp t_ann t_nonl stim msel;
fmdl = rmfield(fmdl,f_to_rem);
imdl.fwd_model = fmdl; imdl_nonl.fwd_model = fmdl;
[imdl_pp(:).fwd_model] = deal(fmdl); [imdl_ann(:).fwd_model] = deal(fmdl);
if distort.large_dist && do_simulation;
    fmdl_fwd = arrayfun(@(x) rmfield(x,f_to_rem), fmdl_fwd);
end

% Less images, keep fmdl so we can easily display images
if distort.large_dist && do_simulation
    c = mat2cell(fmdl_fwd,ones(n_imgs_test,1),1);
    [inh_t(:).fwd_model] = deal(c{:});
else
    [inh_t(:).fwd_model] = deal(rmfield(fmdl_fwd,'coarse2fine'));
end
% Smooth Output of ANN
img_pp = smoothFEMdata( img_pp, 1);
img_ann = smoothFEMdata( img_ann, 1);

%% 6 Display result
txt_titles = { 'Original','One-step GN','PDIPM','ANN','One-step GN + ANN'};
if ~do_simulation; txt_titles = txt_titles(2:end); end
for k=1:1:n_imgs_test
    if do_simulation
        to_disp = [inh_t(k) img_lin(k) img_nonl(k) img_ann(k) img_pp(k)];
    else
        to_disp = [img_lin(k) img_nonl(k) img_ann(k) img_pp(k)];
    end
    [~]=show_multi_fem(to_disp, 'abscissa',txt_titles);
end

%% 7 Display slices
figure; show_slices(img_ann);
figure; show_slices(img_pp)
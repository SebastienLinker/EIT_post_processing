% Post processing 2D EIT bounded problems
% This file intents to train and test ANN for post processing real phantom
% experiments, adjacent current pattern
% Here we do not make use of function newrb or newrbe, known for their poor
% generalization capabilities
%
% (C) 2015/05/22 Sebastien Martin

clearvars; % close all;

%% Parameters
do_simulation = true;
% Tank
n_elecs = 16; % # of electrodes
elec_diam = 0.4; % diameter of one electrode (mm)
bk = 0.0001;

% Targets
n_imgs = 1000;
n_imgs_test = 10;
max_sz= [0.3];
min_sz= [0.15];
min_pos = 0; max_pos = 0.8;% Distance from center
cond_r = [eps 1]; rot_angle = [0 120]; pos_angle = [0 360];
fcn_mk_targ = @createEllipseGeneral; % For training ANN
n_targs = 1;

% Stimulation
opts_stim = {'no_meas_current','no_rotate_meas','do_redundant'};
amp = 0.1; %mA
noisy = false;

% Mesh distortion
distort = true;
distort_tr = false;
dist_type = 'complex'; %'complex'; % 'ellipse';
dist_params.X = [0.5 1]; dist_params.Y = [1 1.5];
dist_params.conformal = [-0.4 0.4];
dist_params.dual_conf = {[-0.32 0],[0 0.32]};

% Inverse solver
prior = @prior_noser;
do_parallel = false;

% ANN
createNN = true;
ann_param.fcn = 'traingdx';
ann_param.goal = 0.0001;
ann_param.min_grad = 1e-10;
ann_param.epochs = 10000;
ann_param.trans = {'radbas','tansig'};
ann_param.max_fail = 500; % # of consecutive iterations without improved performance
ann_param.hid_lay = 500; %hidden layer
ann_param_nopp = ann_param;
ann_param_nopp.trans = {'radbas','tansig'};
ann_param_nopp.hid_lay = 50;

%% 2 Create tank model and inverse model
maxh = 0.05;
imdl_lin = mk_common_model('e2C',n_elecs);
imdl_lin.fwd_model = mdl_normalize(imdl_lin.fwd_model,0);
imdl_lin.solve = @inv_solve_diff_GN_one_step;
fmdl = imdl_lin.fwd_model;
imdl_lin.RtR_prior = prior;
fmdl.jacobian_bkgnd.value = bk;
n_elements = size(fmdl.elems,1);

% Forward models
if noisy
    fmdl.solve = @fwd_solve_add_noise;
    fmdl.add_noise.filter = false;
    % 	fmdl.add_noise.solve = @fwd_solve_higher_order;
    % 	fmdl.approx_type = 'tri6';
    % 	fmdl.system_mat= @system_mat_higher_order;
    % 	fmdl.jacobian = @jacobian_adjoint_higher_order;
end
% Inverse models
imdl_lin.fwd_model = fmdl;
imdl_nonl = imdl_lin;
imdl_nonl.solve = @inv_solve_diff_pdipm;
imdl_nonl.parameters.term_tolerance= 0;
imdl_nonl.parameters.max_iterations = 5;
imdl_nonl.inv_solve_diff_pdipm.norm_data = 1;
imdl_nonl.R_prior = prior;

if distort || distort_tr
    switch dist_type
        case 'ellipse'
            args_dist = {dist_params.X, dist_params.Y};
        case 'complex'
            args_dist = {dist_params.conformal};
        case 'dual_complex'
            args_dist = dist_params.dual_conf(:);
    end
end

%% Create ANNs
if createNN
    %% 3 Randomly creates images
    [conds_NN, area_NN] = fcn_mk_targ(n_imgs,fmdl,bk, max_sz,min_sz, ...
        min_pos, max_pos, rot_angle,pos_angle, cond_r);
    % For multiple targets in a single image
    % 	[ conds_NN, area_NN ] = createMultiTarg( n_imgs, n_targs, 0, fmdl, bk, @createEllipseGeneral,...
    % 						max_sz,min_sz, min_pos,max_pos, rot_angle,pos_angle, cond_r);
    
    % Get the corresponding training data
    [stim,msel]= mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
    fmdl.stimulation = stim;
    fmdl.meas_select = msel;
    imdl_lin.fwd_model = fmdl;
    % Forward and inverse solvers
    if distort_tr
        fmdl_sd = mk_distortion(fmdl, n_imgs, dist_type, args_dist{:});
        for l = 1:1:n_imgs
            [inh_NN(l),~,~,~,t_vh(1,l),t_vi(1,l)] = solveAndNormalize(conds_NN(:,l), fmdl_sd(l), bk, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    else
        [inh_NN,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_NN, fmdl, bk, do_parallel);
    end
    [~,~,t_sol,~,~,~] = solveAndNormalize(0,imdl_lin,bk,do_parallel,[t_vh;t_vi]');
    
    %% 4 Actually creates 1 ANN for PP + 1 ANN with the classic method
    imdl_pp = create_ANN_post_proc(imdl_lin, t_sol, inh_NN, ann_param, bk, 'no_norm', 'do_nodes');
    imdl_ann = create_ANN_inv_solve(imdl_lin, [t_vh;t_vi]', inh_NN, ann_param_nopp, 'do_nodes');
    clear t_vh t_vi t_sol stim msel;
    disp('Trained ANNs');
    
    %% Clean up structures (and save memory)
    f_to_rem = {'stimulation','meas_select'}; %fields to remove
    fmdl = rmfield(fmdl,f_to_rem);
    imdl_lin.fwd_model = rmfield(imdl_lin.fwd_model, f_to_rem);
    imdl_pp.fwd_model = rmfield(imdl_pp.fwd_model, f_to_rem);
    imdl_ann.fwd_model = rmfield(imdl_ann.fwd_model, f_to_rem);
    [inh_NN(:).fwd_model] = deal([]);
else
    disp('Load some pretrained models, please update paths')
    f_to_rem = {'stimulation','meas_select'}; %fields to remove
    load('saved_models.mat','imdl_ann');
    load('saved_models.mat','imdl_pp');
end

%% 5 Randomly generate new images
[conds_test, area_test] = fcn_mk_targ(n_imgs_test,fmdl,bk, max_sz,min_sz, ...
    min_pos, max_pos, rot_angle,pos_angle, cond_r);
xyzr = [area_test(:,4) area_test(:,5) zeros(n_imgs_test,1) area_test(:,1)];
% [conds_test, area_test] = createMultiTarg( n_imgs_test, 2, overlap, fmdl, bk, fcn_targ_test,...
% 						max_sz,min_sz, min_pos,max_pos, rot_angle,pos_angle, cond_r);
% area_test = cell2mat( cellfun(@(x) x(:,1), area_test, 'UniformOutput',false) )';
%% Shape deformation
if distort
    fmdl_sd = mk_distortion(fmdl, n_imgs_test, dist_type, args_dist{:});
else
    fmdl_sd = repmat(fmdl,n_imgs_test,1);
end

%% 6 Apply ANN for post processing and linear solver (comparison purpose)
clear vh_test vi_test img_lin img_nonl img_ann img_pp;
[stim, msel] = mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
if distort && do_simulation
    [fmdl_sd(:).stimulation] = deal(stim);
    [fmdl_sd(:).meas_select] = deal(msel);
end
fmdl.stimulation=stim; fmdl.meas_select=msel;
imdl_lin.fwd_model=fmdl; imdl_pp.fwd_model=fmdl;
imdl_ann.fwd_model = fmdl; imdl_nonl.fwd_model = fmdl;
imdl_nonl.RtR_prior = @prior_noser;
try imdl_nonl = rmfield(imdl_nonl,'R_prior'); catch; end
% Get voltages first, in case we add a random noise, make sure the voltages are always the same
if do_simulation
    if distort
        for l = 1:1:n_imgs_test
            [inh_t(l),~,~,~,t_vh(l),t_vi(l)] = solveAndNormalize(conds_test(:,l), fmdl_sd(l), bk, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    else
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl, bk, do_parallel);
    end
else
    disp('Load the measurement data manually')
    t_vh = vh_ph(k,:); t_vi = vi_ph(k,:);
    inh_t = mk_image(fmdl,conds_test(:,1));
    inh_t.current_params = []; inh_t.info = [];
    inh_t = repmat(inh_t,1,n_imgs_test);
end
% Linear solver, Post-processing, classic ANN, and nonlinear iterative
% methods
[~,~,t_lin,~,~,~] = solveAndNormalize(0,imdl_lin,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_pp,~,~,~] = solveAndNormalize(0,imdl_pp,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_ann,~,~,~] = solveAndNormalize(0,imdl_ann,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_nonl,~,~,~] = solveAndNormalize(0, imdl_nonl, bk, do_parallel,[t_vh;t_vi]');
% Free space
vh_test(:)=t_vh; vi_test(:)=t_vi; img_lin(:)=t_lin;
img_pp(:)=t_pp; img_ann(:)=t_ann;
img_nonl(:)=t_nonl;
clear t_vh t_vi t_lin t_pp t_ann t_nonl stim msel;
fmdl = rmfield(fmdl,f_to_rem); imdl_lin.fwd_model = fmdl;
imdl_pp.fwd_model = fmdl;
imdl_ann.fwd_model = fmdl;
imdl_nonl.fwd_model = fmdl;
if distort && do_simulation
    fmdl_sd = arrayfun(@(x) rmfield(x,f_to_rem), fmdl_sd);
end

% Less images, keep fmdl so we can easily display images
% [img_lin(:).fwd_model] = deal(fmdl); [img_pp(:).fwd_model] = deal(fmdl);
% [img_ann(:).fwd_model] = deal(fmdl); [img_nonl(:).fwd_model] = deal(fmdl);
if distort && do_simulation
    c = mat2cell(fmdl_sd,ones(n_imgs_test,1),1);
    [inh_t(:).fwd_model] = deal(c{:});
else
    [inh_t(:).fwd_model] = deal(fmdl);
end
% Smooth Output of ANN
img_pp = smoothFEMdata( img_pp, 1);
img_ann = smoothFEMdata( img_ann, 1);

%% 6 Display result
for k=1:1:n_imgs_test
    to_disp = [inh_t(k) img_lin(k) img_nonl(k) img_ann(k) img_pp(k)];
    [~]=show_multi_fem(to_disp,'abscissa',{'Original','One step GN','PDIPM','ANN','One-step GN + ANN'}); %,...
    % 					'add_draw',{@drawEllipse,{area_test{1}(1,[4 5 1 2]),'g','LineWidth',1.5}}); %,'force_dir','higher');
end

%% 7 Errors
clear RMS* AR PE RES SD RNG dRES;
img_orig = inh_t;
img_orig_n = normalize_mdl_dir(img_orig,'higher');
img_lin_n = normalize_mdl_dir(img_lin,'higher');
img_ann_n = normalize_mdl_dir(img_ann,'higher');
img_pp_n = normalize_mdl_dir(img_pp,'higher');
img_nonl_n = normalize_mdl_dir(img_nonl,'higher');
xyzr = repmat([0;0;0;0],1,n_imgs_test);
if ~do_simulation
    if iscell(area_test)
        area_t = area_test{1}(:,1)'; %Used for weighted RMS only, it's fucked with two targets
    end
    area_t = repmat(area_t,n_imgs_test,1);
else
    area_t = area_test;
end
imgs_GREIT = [img_orig; img_lin; img_nonl; img_ann; img_pp];
for k = 1:1:n_imgs_test
    ell_pos = [area_t(k,[1 2]) 0 area_t(k,[4 5])];
    RMS_lin(k) = calc_errorRMS(img_orig_n(k), img_lin_n(:,k), 'area');
    RMS_ann(k) = calc_errorRMS(img_orig_n(k), img_ann_n(:,k), 'area');
    RMS_pp(k) = calc_errorRMS(img_orig_n(k), img_pp_n(:,k), 'area');
    RMS_nonl(k) = calc_errorRMS(img_orig_n(k), img_nonl_n(:,k), 'area');
    [~, AR(:,k), PE(:,k), dRES(:,k), SD(:,k), RNG(:,k)] = cmp_multi_GREIT( imgs_GREIT(2:end,k), img_orig(k), xyzr(:,k));
end
RMS_errors = [RMS_lin; RMS_nonl; RMS_ann; RMS_pp];

disp('One column per test image')
disp('For each image, display errors with: Linear method, NonLinear, ANN, and Post processing')
disp('RMS Errors:')
disp(RMS_errors)
disp('')
disp('AR')
disp(AR)
disp('')
disp('PE')
disp(PE)
disp('')
disp('dRES')
disp(dRES)
disp('')
disp('SD')
disp(SD)
disp('')
disp('RNG')
disp(RNG)

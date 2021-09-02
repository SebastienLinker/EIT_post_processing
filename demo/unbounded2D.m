% Post processing 2D EIT unbounded problems
% This file intents to train and test ANN for post processing real phantom
% experiments, adjacent current pattern
%
% (C) 2015/05/22 Sebastien Martin

clearvars; % close all;

%% Parameters
do_simulation = true;
% Tank
n_elecs = 8; % # of electrodes
elec_diam = 0.2; % diameter of one electrode (mm)
bk = 0.0001;
bk_est = 0.0001;

probe_rad = 1;
sz_small_mdl = 15;

% Targets
n_imgs = 1000;
n_imgs_test = 10;
max_sz= [4 7];
min_sz= [2 5];
min_pos = 0; max_pos = 3;% Distance from center
cond_r = [0.9 100]; rot_angle = [-180 180]; pos_angle = [0 360];
pos_ang_test = 0;
fcn_mk_targ = @createEllipseGeneral; % For training ANN
n_targs = 1;
fcn_targ_test = @createEllipseGeneral; % Testing
overlap = false;

% Stimulation
opts_stim = {'no_meas_current','no_rotate_meas','do_redundant'};
amp = 0.1; %mA
noisy = false;

% Mesh distortion
distort = false;
distort_tr = false;
dist_type = 'ellipse'; %'complex'; % 'ellipse';
dist_params.X = [0.5 1]; dist_params.Y = [1 1.5];
dist_params.conformal = [-0.4 0.4];
dist_params.dual_conf = {[-0.32 0],[0 0.32]};

% Inverse solver
prior = @prior_laplace;
p_type = 'R_prior';
do_parallel = false;

% ANN
createNN = true; % True if you want to use a NN, complete
ann_param.fcn = 'traingdx';
ann_param.goal = 0.0001;
ann_param.min_grad = 1e-10;
ann_param.epochs = 1000;
ann_param.trans = {'radbas','tansig'};
ann_param.max_fail = 500; % # of consecutive iterations without improved performance
ann_param.hid_lay = 1000; %hidden layer
ann_param_nopp = ann_param;
ann_param_nopp.trans = {'radbas','tansig'};
ann_param_nopp.hid_lay = 800;

%% 2 Create tank model and inverse model
maxh = 0.05;
imdl_lin = rmfield( mk_common_model('e2C',n_elecs), 'RtR_prior');
fmdl = mk_unbounded2D_phantom(.75, n_elecs);
imdl_lin.fwd_model = fmdl;
imdl_lin.(p_type) = prior;
fmdl.jacobian_bkgnd.value = bk;
n_elements = size(fmdl.elems,1);
% Forward models
if noisy
    fmdl.solve = @fwd_solve_add_noise;
    fmdl.add_noise.filter = true;
    fmdl.add_noise.solve = @fwd_solve_higher_order;
    fmdl.approx_type = 'tri6';
    fmdl.system_mat= @system_mat_higher_order;
    fmdl.jacobian = @jacobian_adjoint_higher_order;
end
% Inverse models
imdl_lin.fwd_model = fmdl;
imdl_nonl = imdl_lin;
imdl_nonl.solve = @inv_solve_diff_pdipm;
imdl_nonl.parameters.term_tolerance= 0;
imdl_nonl.parameters.max_iterations = 5;
imdl_nonl.inv_solve_diff_pdipm.norm_data = 1;
imdl_nonl = rmfield(imdl_nonl,p_type);
imdl_nonl.RtR_prior = @prior_noser;

if distort
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
    %% 3.5 Get the corresponding training data
    [stim,msel]= mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
    fmdl.stimulation = stim; fmdl.meas_select = msel; imdl_lin.fwd_model = fmdl;
    % Forward and inverse solvers
    if distort_tr
        fmdl_sd = mk_distortion(fmdl, n_imgs, dist_type, args_dist{:});
        for l = 1:1:n_imgs
            [inh_NN(l),~,~,~,t_vh(l),t_vi(l)] = solveAndNormalize(conds_NN(:,l), fmdl_sd(l), bk, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    else
        [inh_NN,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_NN, fmdl, bk, do_parallel);
    end
    tmp = t_vh; t_vh = t_vi; t_vi = tmp; clear tmp;
    [~,~,t_sol,~,~,~] = solveAndNormalize(0,imdl_lin,bk_est,do_parallel,[t_vh;t_vi]');
    
    %% 4 Actually creates 1 ANN for PP + 1 ANN with the classic method
    imdl_pp = create_ANN_post_proc(imdl_lin, t_sol, inh_NN, ann_param, bk_est, 'no_norm', 'do_nodes');
    imdl_ann = create_ANN_inv_solve(imdl_lin, [t_vh;t_vi]', inh_NN, ann_param_nopp, 'do_nodes');
    vh_NN=t_vh; vi_NN=t_vi; sol_NN=t_sol; clear t_vh t_vi t_sol stim msel;
    disp(['Trained ANNs']);
    
    %% Clean up structures
    f_to_rem = {'stimulation','meas_select'}; %fields to remove
    fmdl = rmfield(fmdl,f_to_rem);
    imdl_lin.fwd_model = rmfield(imdl_lin.fwd_model, f_to_rem);
    imdl_pp.fwd_model = rmfield(imdl_pp(1).fwd_model, f_to_rem);
    imdl_ann.fwd_model = rmfield(imdl_ann(1).fwd_model, f_to_rem);
    [inh_NN(:).fwd_model] = deal([]);
    [sol_NN(:).fwd_model] = deal([]);
    clear conds_NN inh_NN_n sol_NN_n sol_NN; % Keep voltages in case of noise
else
    % Load data
end

%% 5 Randomly generate new images
pos_targ = [3.0:0.1:3.9]; % pos_targ = [4.1:0.1:8.0];
clear conds_test area_test;
for k = 1:1:n_imgs_test
    % Make images one by one, cause we have different forward models
    [conds_test(:,k), area_test(k,:)] = fcn_mk_targ(1,fmdl,bk, [3 6], [3 6], ...
        pos_targ(k), pos_targ(k), [0 0],pos_ang_test, 0.9);
end
xyzr = [area_test(:,4) area_test(:,5) zeros(n_imgs_test,1) area_test(:,1)];

%% Shape deformation
if distort
    fmdl_sd = mk_distortion(fmdl, n_imgs_test, dist_type, args_dist{:});
else
    fmdl_sd = [];
end

%% 5.5 If phantom measurements, load data
if ~do_simulation
    disp('Please load your own data')
    
    % Simulate initial target
    n_imgs_test = n_frames;
    img_sz = {[2 2]};
    dist = {1+2+7.5/2.5}; ang = {0}; cond_targ = {0.9};
    [conds_test, area_test] = createEllipseGeneral(1, fmdl, bk,img_sz{1},img_sz{1},dist{1},dist{1},[0 0],ang{1},cond_targ{1});
    inh_t = repmat( mk_image(fmdl, conds_test), 1, n_imgs_test);
    area_test = repmat(area_test, n_imgs_test,1);
end

%% 6 Apply ANN for post processing and linear solver (comparison purpose)
clear vh_test vi_test img_lin img_nonl img_ann img_pp;
[stim, msel] = mk_stim_patterns(n_elecs,1,'{ad}','{ad}',opts_stim,amp);
if distort && do_simulation
    [fmdl_sd(:).stimulation] = deal(stim);
    [fmdl_sd(:).meas_select] = deal(msel);
end
fmdl.stimulation=stim;
fmdl.meas_select=msel;
imdl_lin.fwd_model=fmdl;
imdl_pp.fwd_model=fmdl;
imdl_ann.fwd_model = fmdl;
imdl_nonl.fwd_model = fmdl;
% Get voltages first, in case we add a random noise, make sure the voltages are always the same
if do_simulation
    c_t = conds_test;
    bk_n = 0.2;
    c_t(c_t==bk) = bk_n;
    if distort
        for l = 1:1:n_imgs_test
            [inh_t(l),~,~,~,t_vh(l),t_vi(l)] = solveAndNormalize(conds_test(:,l), fmdl_sd(l), bk, do_parallel);
            [~,~,~,~,t_vhN(l),t_viN(l)] = solveAndNormalize(c_t(:,l), fmdl_sd(l), bk_n, do_parallel);
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
        t_vhN = t_vhN(randperm(length(t_vhN)));
    else
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl, bk, do_parallel);
        [~,~,~,~,t_vhN,t_viN] = solveAndNormalize(c_t, fmdl, bk_n, do_parallel);
    end
else
    disp('Load measurement data')
end
% Linear solver, Post-processing, classic ANN
[~,~,t_lin,~,~,~] = solveAndNormalize(0,imdl_lin,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_ann,~,~,~] = solveAndNormalize(0,imdl_ann,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_pp,~,~,~] = solveAndNormalize(0,imdl_pp,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_nonl,~,~,~] = solveAndNormalize(0, imdl_nonl, bk_est, do_parallel,[t_vhN;t_viN]');
vh_test=t_vh;
vi_test=t_vi;
img_lin=t_lin;
img_pp=t_pp;
img_ann=t_ann;
img_nonl=t_nonl;
clear t_vh t_vi t_vhN t_viN c_t t_lin t_pp t_ann t_nonl stim msel;
fmdl = rmfield(fmdl,f_to_rem);
imdl_lin.fwd_model = fmdl;
imdl_pp.fwd_model = fmdl;
imdl_ann.fwd_model = fmdl;
imdl_nonl.fwd_model = fmdl;
if distort && do_simulation
    fmdl_sd = arrayfun(@(x) rmfield(x,f_to_rem), fmdl_sd);
end
% Less images, keep fmdl so we can easily display images
if distort && do_simulation
    c = mat2cell(fmdl_sd,ones(n_imgs_test,1),1);
    [inh_t(:).fwd_model] = deal(c{:});
else
    [inh_t(:).fwd_model] = deal(fmdl);
end
% Smooth Output of ANN
img_pp = smoothFEMdata( img_pp, 1);
img_ann = smoothFEMdata( img_ann, 1);

%% 7 Display result
for k=1:1:n_imgs_test
    to_disp = [repmat(inh_t(1,k),1,1) ...
        img_lin(k) img_nonl(k) img_ann(k) img_pp(k)];
    [~]=show_multi_fem(to_disp,'abscissa',{'Original','One step GN','PDIPM','ANN','ANN PP'},...
        'add_draw',{@drawEllipse,{area_test(k,[4 5 1 2]),'g','LineWidth',1.5}}); %,'force_dir','higher');
    clear to_disp;
end

%% 8 'Clear' display
for k=1:1:n_imgs_test
    figure;
    subtightplot(2,2,1);
    show_fem(img_lin(k)); axis off; axis(sz_small_mdl.*[-1 1 -1 1]); title('One step GN');
    hold on; drawEllipse(area_test(k,[4 5 1 2]),'g','LineWidth',1.5);
    subtightplot(2,2,2);
    show_fem(img_nonl(k)); axis off; axis(sz_small_mdl.*[-1 1 -1 1]); title('PDIPM');
    hold on; drawEllipse(area_test(k,[4 5 1 2]),'g','LineWidth',1.5);
    subtightplot(2,2,3);
    show_fem(img_pp(k)); axis off; axis(sz_small_mdl.*[-1 1 -1 1]); title('ANN');
    hold on; drawEllipse(area_test(k,[4 5 1 2]),'g','LineWidth',1.5);
    subtightplot(2,2,4);
    show_fem(img_ann(k)); axis off; axis(sz_small_mdl.*[-1 1 -1 1]); title('One step GN + ANN');
    hold on; drawEllipse(area_test(k,[4 5 1 2]),'g','LineWidth',1.5);
end

%% 9 RMS, GREIT Errors
clear RMS* AR PE RES SD RNG dRES n_err err_dist err_area_pct th;
inh_t_n = normalize_mdl_dir(inh_t,'higher');

% Estimate background
img_lin_n = normalize_mdl(img_lin, 'graphics','higher');
img_nonl_n = normalize_mdl(img_nonl, 'graphics','higher');
img_ann_n = normalize_mdl(img_ann, 'graphics','higher');
img_pp_n = normalize_mdl(img_pp, 'graphics','higher');
xyzr = [area_test(:,[4 5]), zeros(size(area_test,1),1), area_test(:,1)]';
if ~do_simulation
    xyzr = repmat(xyzr,1,n_imgs_test);
    area_test = repmat(area_test,n_imgs_test,1);
end

imgs_err = [img_lin; img_nonl; img_ann; img_pp];
imgs_err_n = [img_lin_n; img_nonl_n; img_ann_n; img_pp_n];
imgs_err_greit = [img_lin; img_nonl_n; img_ann_n; img_pp_n];

% GREIT
bnd_gr = [zeros(1,n_imgs_test); xyzr(1,:)+4.5; -9.*ones(1,n_imgs_test); 9.*ones(1,n_imgs_test)];
bnd_gr = [xyzr(1,:)-6; xyzr(1,:)+4.5; -9.*ones(1,n_imgs_test); 9.*ones(1,n_imgs_test)];
greit_func = @(x) cmp_multi_GREIT( imgs_err_greit(:,x), inh_t(x), xyzr(:,x), 'no_norm',32,bnd_gr(:,x));

% Distance error
n_err = zeros(size(imgs_err')); err_dist = zeros(size(imgs_err'));
err_area_pct = zeros(size(imgs_err')); th = zeros(size(imgs_err'));

for k = 1:1:n_imgs_test
    ell_pos = [area_test(k,[1 2]) 0 area_test(k,[4 5])];
    % RMS
    RMS_err(:,k) = calc_errorRMS(inh_t_n(k), imgs_err_n(:,k), 'area');
    RMS_err_w(:,k) = calc_errorRMS(inh_t_n(k), imgs_err_n(:,k), 'area', ...
        'half_ellipse',[area_test(k,1:2)*2, area_test(k,3:5)], 'weight',2/3);
    RMS_err_ell(:,k) = calc_errorRMS(inh_t_n(k), imgs_err_n(:,k), 'area', ...
        'ellipse',[area_test(k,1:2)*2, area_test(k,3:5)], 'weight',2/3);
    % GREIT
    [~, AR(:,k), PE(:,k), dRES(:,k), SD(:,k), RNG(:,k)] = greit_func(k);
    % Distance error
    [~, n_err(k,:), err_dist(k,:), err_area_pct(k,:), th(k,:)] = ...
        area_error_unbounded( inh_t_n(k), imgs_err_n(:,k), area_test(k,:));
end
clear inh_t_n img_lin_n img_nonl_n img_ann_n img_pp_n imgs_err_n imgs_err;

disp('Display order (by row): Linear solver, PDIPM, ANN, post processing')
disp('Target position')
disp(pos_targ)
fprintf('RMS error, focus on ellipse:\n')
disp(RMS_err_ell)
fprintf('AR:\n')
disp(AR)
fprintf('Position error:\n')
disp(PE)
fprintf('Difference of resolution:\n')
disp(dRES)
fprintf('Shape deformation:\n')
disp(SD)
fprintf('Ringing:\n')
disp(RNG)
fprintf('Area error in distance:\n')
disp(err_dist')
fprintf('Area error in distance (percentage):\n')
disp(err_area_pct')

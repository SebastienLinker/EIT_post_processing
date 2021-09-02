% Post processing 3D EIT unbounded problems
% This file trains and tests ANN for post processing real phantom
% experiments, adjacent current pattern
%
% (C) 2015/11/19 Sebastien Martin

clearvars; % close all;

%% Parameters
do_simulation = true;
% Tank
n_elecs = 8; % # of electrodes
n_layers = 4;
bk = 0.0001; 1.5e-4; 1;
bk_est = 1;
c_targ = 1; 2.5e-4; 2;
probe_rad = 1;
sz_small_mdl = 15;

% Targets
n_imgs = 1000;
n_imgs_test = 10;
max_sz= [4 7 10];
min_sz= [2 4 8];
min_pos = 0; max_pos = 4; %[8 -5 1];% Distance from center
cond_r = [eps 1]; rot_angle = [0 180]; pos_angle = [-180 180];
pos_ang_test = 0;
fcn_mk_targ = @createEllipseGeneral; % For training ANN
% fcn_mk_targ = @createCylinder;
% fcn_mk_targ = @createRandom;
n_targs = 1; fcn_targ_test = @createEllipseGeneral; overlap = false; % Testing
% fcn_targ_test = @createCylinder;

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
do_parallel = false;

% ANN
createNN = true; % True if you want to use a NN, complete
ann_param.fcn = 'traingdx';
ann_param.goal = 0.0001;
ann_param.min_grad = 1e-10;
ann_param.epochs = 5000;
ann_param.trans = {'radbas','tansig'};
ann_param.max_fail = 500; % # of consecutive iterations without improved performance
ann_param.hid_lay = 1000; %hidden layer
ann_param_nopp = ann_param;
ann_param_nopp.trans = {'radbas','tansig'};
ann_param_nopp.hid_lay = 1000;

%% 2 Create tank model and inverse model
maxh = 0.05;
% 	fmdl(k) = createBounded2Dtank_1cyl_obj( n_elec, maxh, [X(k,:), Y(k,:), sz_targ(k,:)] );
imdl_lin = rmfield( mk_common_model('e2C',n_elecs), 'RtR_prior');
imdl_lin.solve = @inv_solve_diff_GN_one_step;
% [fmdl, fmdl_s] = createUnbounded3Dphantom(n_layers, sz_small_mdl);
[fmdl, fmdl_s] = createUnbounded3D_simple(n_layers, sz_small_mdl);
imdl_lin.fwd_model = fmdl;
imdl_lin.R_prior = prior;
fmdl.jacobian_bkgnd.value = bk;
n_elements = size(fmdl.elems,1);
% Forward models
if noisy
    fmdl.solve = @fwd_solve_add_noise;
    fmdl.add_noise.filter = true;
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
imdl_nonl.R_prior = prior; @prior_noser;


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
    [stim,msel]= mk_stim_patterns(n_elecs,n_layers,'{ad}','{ad}',opts_stim,amp);
    fmdl.stimulation = stim; fmdl.meas_select = msel; imdl_lin.fwd_model = fmdl;
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
    [~,~,t_sol,~,~,~] = solveAndNormalize(0,imdl_lin,bk_est,do_parallel,[t_vh;t_vi]');
    
    %% 4 Actually creates 1 ANN for PP + 1 ANN with the classic method
    [t_sol(:).fwd_model] = deal(fmdl_s);
    [inh_NN(:).fwd_model] = deal(fmdl_s);
    imdl_pp = create_ANN_post_proc(imdl_lin, t_sol, inh_NN, ann_param, bk_est, 'no_norm', 'do_nodes', 'do_c2f');
    imdl_lin_s = setfield(imdl_lin,'fwd_model',fmdl_s);
    imdl_ann = create_ANN_inv_solve(imdl_lin_s, [t_vh;t_vi]', inh_NN, ann_param_nopp, 'do_nodes');
    clear t_vh t_vi t_sol stim msel imdl_lin_s;
    disp(['Trained ANNs']);
    
    %% Clean up structures (and save lots of memory, especially in case of CS-EIT)
    f_to_rem = {'stimulation','meas_select'}; %fields to remove
    fmdl = rmfield(fmdl,f_to_rem);
    imdl_lin.fwd_model = rmfield(imdl_lin.fwd_model, f_to_rem);
    imdl_pp.fwd_model = rmfield(imdl_pp.fwd_model, f_to_rem);
    clear conds_NN inh_NN_n sol_NN_n vh_NN vi_NN;
else
    % LOAD YOUR PRE TRAINED MODELS HERE
end

%% 5 Randomly generate new images
pos_targ = [4.1:0.1:8.0];
clear conds_test area_test inh_t_s;
for k = 1:1:n_imgs_test
    % Make images one by one, cause we have different forward models
    [conds_test(:,k), area_test(k,:)] = fcn_mk_targ(1,fmdl,bk, [3 6 9], [3 6 9], ...
        pos_targ(k), [pos_targ(k) 0 0], [0 0],pos_ang_test, 0.9);
    tmp = fcn_mk_targ(1,fmdl_s,bk, [3 6 9],[3 6 9], pos_targ(k),[pos_targ(k) 0], [0 0],pos_ang_test, 0.9);
    tmp = setfield( mk_image(fmdl_s,tmp), 'current_params',[]);
    inh_t_s(k) = orderfields( setfield(tmp,'info',struct('error',NaN)) ); clear tmp;
end
xyzr = [area_test(:,4) area_test(:,5) zeros(n_imgs_test,1) area_test(:,1)];

%% Shape deformation
if distort
    fmdl_sd = mk_distortion(fmdl, n_imgs_test, dist_type, args_dist{:});
else
end

%% 5.5 If phantom measurements, load data
if ~do_simulation
    disp('Please load your measurement data')
end

%% 6 Apply ANN for post processing and linear solver (comparison purpose)
clear vh_test vi_test img_lin img_nonl img_ann img_pp;
[stim,msel]= mk_stim_patterns(n_elecs,n_layers,'{ad}','{ad}',opts_stim,amp);
if distort && do_simulation
    [fmdl_sd(:).stimulation] = deal(stim); [fmdl_sd(:).meas_select] = deal(msel);
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
        end
        t_vh = t_vh(randperm(length(t_vh))); % Mix-up the homogeneous data
    else
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl, bk, do_parallel);
    end
else
    t_vh = vh_ph; t_vi = vi_ph;
    inh_t = mk_image(fmdl,conds_test(:,1));
    inh_t.current_params = [];
    inh_t.info = [];
    inh_t = repmat(inh_t,1,n_imgs_test);
end
% Linear solver, Post-processing, classic ANN
[~,~,t_lin,~,~,~] = solveAndNormalize(0,imdl_lin,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_ann,~,~,~] = solveAndNormalize(0,imdl_ann,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_pp,~,~,~] = solveAndNormalize(0,imdl_pp,bk_est,do_parallel,[t_vh;t_vi]');
[~,~,t_nonl,~,~,~] = solveAndNormalize(0, imdl_nonl, bk_est, do_parallel,[t_vh;t_vi]');
% Free space
img_lin = t_lin;
img_pp=t_pp;
img_ann=t_ann;
img_nonl=t_nonl;
clear t_vh t_vi t_vhN t_viN t_lin t_pp t_ann t_nonl stim msel;
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
    inh_t = inh_t(1,:); % They are all the same
else
    [inh_t(:).fwd_model] = deal(fmdl);
end

% Minimize the size
[img_pp(:).fwd_model] = deal(fmdl_s);
[img_ann(:).fwd_model] = deal(fmdl_s);
[img_lin(:).fwd_model] = deal(fmdl_s);
[img_nonl(:).fwd_model] = deal(fmdl_s);

% Smooth Output of ANN
img_pp = smoothFEMdata( img_pp, 1);
img_ann = smoothFEMdata( img_ann, 1);

%% 7 Display result
for k=1:1:n_imgs_test
    to_disp = [repmat(inh_t_s(1,k),1) img_lin(k) img_nonl(k) img_ann(k) img_pp(k)];
    % 	[~]=show_multi_fem(to_disp,'abscissa',{'Original','One step GN','PDIPM','ANN','ANN PP'},...
    % 					'add_draw',{@draw_ellipsoid_in_fem,{area_test(k,[6 7 8 1 2 3 4 5 5])}});
    clear to_disp;
end

%% 8 'Clear' display
for k=1:1:n_imgs_test
    figure;
    subtightplot(2,2,1);
    show_fem(img_lin(k));
    axis off;
    title('One step GN');
    subtightplot(2,2,2);
    show_fem(img_nonl(k));
    axis off;
    title('PDIPM');
    subtightplot(2,2,3);
    show_fem(img_ann(k));
    axis off;
    title('ANN');
    subtightplot(2,2,4);
    show_fem(img_pp(k));
    axis off;
    title('One step GN + ANN');
end

%% Show slices
slice_Z = 1.5; n_pts = 128;
for k= 1:1:n_imgs_test
    pos_targ = [area_test(k,6); -area_test(k,7)] / sz_small_mdl*(n_pts/2)+(n_pts/2);
    sz_targ = area_test(k,[1 2])'/sz_small_mdl*(n_pts/2);
    to_disp = [img_lin(k), img_nonl(k), img_ann(k), img_pp(k)];
    [to_disp(:).calc_colours] = deal( struct('npoints',n_pts) );
    [~]=show_multi_fem(to_disp,'abscissa',{'One-step GN','PDIPM','ANN','One step GN + ANN'},...
        'slice',[Inf Inf slice_Z],...
        'add_draw',{@drawEllipse,{[pos_targ; sz_targ]','g','LineWidth',1.5}});
    clear to_disp;
end

%% 9 Mapping into 2D model (test)
img_lin_n = rmfield( normalize_mdl_dir(img_lin,'higher'), 'calc_colours');
img_nonl_n = normalize_mdl(img_nonl, 'graphics','higher');
img_ann_n = normalize_mdl(img_ann, 'graphics','higher');
img_pp_n = normalize_mdl(img_pp, 'graphics','higher');
inh_t_n = normalize_mdl_dir(inh_t_s,'higher');

fmdl2D = mk_unbounded2D_phantom(1, n_elecs);
fmdl2D.mk_coarse_fine_mapping.z_depth = 0.2;
fmdl2D.mk_coarse_fine_mapping.f2c_offset = [0 0 1.5];
fmdl2D.coarse2fine = mk_coarse_fine_mapping( fmdl_s, fmdl2D )';
img_ann_2D = img_ann_n;
[img_ann_2D(:).fwd_model] = deal(fmdl2D);
img_pp_2D = img_pp_n;
[img_pp_2D(:).fwd_model] = deal(fmdl2D);
img_lin_2D = img_lin_n;
[img_lin_2D(:).fwd_model] = deal(fmdl2D);
img_nonl_2D = img_nonl_n;
[img_nonl_2D(:).fwd_model] = deal(fmdl2D);
inh_t_2D = inh_t_n;
[inh_t_2D(:).fwd_model] = deal(fmdl2D);
for k=1:1:n_imgs_test
    to_disp = [img_lin_2D(k), img_nonl_2D(k), img_ann_2D(k), img_pp_2D(k)];
    % 	[~]=show_multi_fem(to_disp, 'axis', [-15 15 -15 15], ...
    % 				'ordinate',{'One-step GN','PDIPM','ANN','One step GN + ANN'});
    clear to_disp;
end
imgs_err_2D_n = [img_lin_2D; img_nonl_2D; img_ann_2D; img_pp_2D];

%% 10 RMS, GREIT Errors
clear RMS* AR PE RES SD RNG dRES n_err err_dist err_area_pct th;

inh_t_n = normalize_mdl_dir(inh_t_s,'higher');

% Estimate background
img_lin_n = normalize_mdl(img_lin, 'graphics','higher');
img_nonl_n = normalize_mdl(img_nonl, 'graphics','higher');
img_ann_n = normalize_mdl(img_ann, 'graphics','higher');
img_pp_n = normalize_mdl(img_pp, 'graphics','higher');

xyzr = [area_test(:,[6 7]), 1*ones(size(area_test,1),1), area_test(:,1)]';
if ~do_simulation
    xyzr = repmat(xyzr,1.5,n_imgs_test);
    area_test = repmat(area_test,n_imgs_test,1);
end

imgs_err = [img_lin; img_nonl; img_ann; img_pp];
imgs_err_n = [img_lin_n; img_nonl_n; img_ann_n; img_pp_n];

% GREIT
bnd_gr = [xyzr(1,:)-9; xyzr(1,:)+9; -9.*ones(1,n_imgs_test); 9.*ones(1,n_imgs_test)];
% Use this one for 3D
% greit_func = @(x) cmp_multi_3D_GREIT( imgs_err_n(:,x), inh_t_n(x), 'no_norm',32,bnd_gr(:,x));
greit_func = @(x) cmp_multi_3D_GREIT( imgs_err_n(:,x), inh_t_n(x), false, 32, [-2 12 -7 7]);
% 2D GREIT
% greit_func = @(x) cmp_multi_GREIT( imgs_err_n(:,x), inh_t_n(x), xyzr(:,x), 'no_norm',32,bnd_gr(:,x));

% Distance error
n_err = zeros(size(imgs_err'));
err_dist = zeros(size(imgs_err'));
err_area_pct = zeros(size(imgs_err'));
th = zeros(size(imgs_err'));

for k = 1:1:n_imgs_test
    ell_pos = [area_test(k,[1 2]) 0 area_test(k,[4 5])];
    % RMS
    RMS_err(:,k) = calc_errorRMS(inh_t_n(k), imgs_err_n(:,k), 'area');
    RMS_err_ell(:,k) = calc_errorRMS(inh_t_n(k), imgs_err_n(:,k), 'area', ...
        'ellipse',[area_test(k,1:3)*2, area_test(k,4:8)], 'weight',2/3);
    % GREIT
    [~, AR(:,k), PE(:,k), dRES(:,k), SD(:,k), RNG(:,k)] = greit_func(k);
    % Distance error
    [~, n_err(k,:), err_dist(k,:), err_area_pct(k,:), th(k,:)] = ...
        area_error_unbounded( inh_t_2D(k), imgs_err_2D_n(:,k), area_test(k,[1 2 4 6 7 9]));
end

%% Display errors
disp('Display order (by row): Linear solver, PDIPM, ANN, post processing')
disp('Target position (distance from target''s center to model''s center)')
disp(area_test(:,6)')
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

% Post processing 3D EIT bounded problems
% This file intents to train and test ANN for post processing real phantom
% experiments, adjacent current pattern
%
% (C) 2015/05/22 Sebastien Martin

clearvars; close all;

%% Parameters
do_simulation = true;
% Tank
n_elecs = 8; % # of electrodes per layer

n_layers = 4;
elec_diam = 0.4; % diameter of one electrode (mm)
bk = 0.0001;
calcRM = 'getRM'; % 'getRM' or 'noRM'

% Targets
n_imgs = 100;
n_imgs_test = 10;
max_sz= [0.3];
min_sz= [0.1]; 
min_pos = 0; max_pos = 0.8;% Distance from center
cond_r = [0.2 1]; rot_angle = [0 180]; pos_angle = [0 360];
fcn_mk_targ = @createCylinder; % For training ANN

n_targs = 2; overlap = false; % Testing

% Stimulation
opts_stim = {'no_meas_current','no_rotate_meas','do_redundant'};
amp = 8; %mA
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
ann_param.epochs = 100000;
ann_param.trans = {'radbas','tansig'};
ann_param.max_fail = 5000; % # of consecutive iterations without improved performance
ann_param.hid_lay = 250; %hidden layer
ann_param_nopp = ann_param;


%% 2 Create tank model and inverse model
maxh = 0.05;
% 	fmdl(k) = createBounded2Dtank_1cyl_obj( n_elec, maxh, [X(k,:), Y(k,:), sz_targ(k,:)] );
imdl_lin = mk_common_model('e2C',16);
fmdl = mk_bounded3D_phantom_cross_sec();
fmdl = mk_bounded3D_simple();
fmdl.jacobian_bkgnd.value = bk;
% [stim, meas_sel] = mk_multi_pattern(n_elecs,1,'{phantom_ad}',[0,1],CS_step,'LHC',options_stim,amp);
% fmdl.stimulation = stim; fmdl.meas_select = meas_sel;
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
	if n_targs==1
		[conds_NN, area_NN] = fcn_mk_targ(n_imgs,fmdl,bk, max_sz,min_sz, ...
								min_pos, max_pos, rot_angle,pos_angle, cond_r);
	else
		[conds_NN, area_NN] = createMultiTarg(n_imgs,n_targs,overlap,fmdl,bk, fcn_mk_targ,...
				max_sz,min_sz,min_pos, max_pos, rot_angle,pos_angle, cond_r);
    end
	%% 3.5 Get the corresponding training data
    % Generate stimulations for this CS step delta alpha
    [stim,msel]= mk_stim_patterns(n_elecs,n_layers,'{ad}',[0,1],opts_stim,amp);

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
    [~,~,t_sol,~,~,~] = solveAndNormalize(0,imdl_lin,bk,do_parallel,[t_vh;t_vi]');

    %% 4 Actually creates 1 ANN for PP + 1 ANN with the classic method
    imdl_pp = create_ANN_post_proc(imdl_lin, t_sol, inh_NN, ann_param, bk, ...
                                        'no_norm','do_nodes',calcRM);
    imdl_ann = create_ANN_inv_solve(imdl_lin, [t_vh;t_vi]', inh_NN, ann_param_nopp, 'do_nodes');
    % Keep voltages for later use?
    vh_NN=t_vh; vi_NN=t_vi; sol_NN=t_sol;
    clear t_vh t_vi t_sol stim msel;
    disp(['Trained ANNs']);
	
	%% Clean up structures (and save lots of memory)
	f_to_rem = {'stimulation','meas_select'}; %fields to remove
	fmdl = rmfield(fmdl,f_to_rem);
	imdl_lin.fwd_model = rmfield(imdl_lin.fwd_model, f_to_rem);
	imdl_pp.fwd_model = rmfield(imdl_pp.fwd_model, f_to_rem);
	imdl_ann.fwd_model = rmfield(imdl_ann(1).fwd_model, f_to_rem);
	[inh_NN(:).fwd_model] = deal([]); 
	[sol_NN(:).fwd_model] = deal([]); 
else
    % Load pre trained models
end

%% 5 Randomly generate new images
if n_targs==1
	[conds_test, area_test] = fcn_mk_targ(n_imgs_test,fmdl,bk, max_sz,min_sz, ...
								min_pos, max_pos, rot_angle,pos_angle, cond_r);
else
	[conds_test, area_test] = createMultiTarg(n_imgs_test,n_targs,overlap,fmdl,bk, fcn_mk_targ,...
					max_sz,min_sz,min_pos, max_pos, rot_angle,pos_angle, cond_r);
% 	area_test = cell2mat(cellfun(@(x) x(:,1),area_test,'UniformOutput',false))'; % no crash
end

xyzr = [area_test(:,4) area_test(:,5) zeros(n_imgs_test,1) area_test(:,1)];
% [conds_test, area_test] = createMultiTarg( n_imgs_test, n_targs, overlap, fmdl, bk, fcn_targ_test,...
% 						max_sz,min_sz, min_pos,max_pos, rot_angle,pos_angle, cond_r);
%% Shape deformation
if distort
	fmdl_sd = mk_distortion(fmdl, n_imgs_test, dist_type, args_dist{:});
	t_vh_order = randperm(n_imgs_test);
else
    fmdl_sd = repmat(fmdl,n_imgs_test,1);
end
							
%% 5.5 If phantom measurements, load data
if ~do_simulation
    disp('Please load measurement data')
end

%% 6 Apply ANN for post processing and linear solver (comparison purpose)
clear vh_test vi_test img_lin img_nonl img_ann img_pp;

[stim, msel] = mk_stim_patterns(n_elecs,n_layers,'{ad}',[0,1],opts_stim,amp);
if distort && do_simulation
    [fmdl_sd(:).stimulation] = deal(stim); 
    [fmdl_sd(:).meas_select] = deal(msel);
end
fmdl.stimulation=stim; fmdl.meas_select=msel;
imdl_lin.fwd_model=fmdl;
imdl_pp.fwd_model=fmdl;
imdl_ann.fwd_model = fmdl;
imdl_nonl.fwd_model = fmdl;
% Get voltages first, in case we add a random noise, make sure the voltages are always the same
if do_simulation
    if distort
        for l = 1:1:n_imgs_test
            [inh_t(l),~,~,~,t_vh(l),t_vi(l)] = solveAndNormalize(conds_test(:,l), fmdl_sd(l), bk, do_parallel);
        end
        t_vh = t_vh(t_vh_order); % Mix-up the homogeneous data
    else
        [inh_t,~,~,~,t_vh,t_vi] = solveAndNormalize(conds_test, fmdl, bk, do_parallel);
    end
else
    t_vh = vh_ph(k,:); t_vi = vi_ph(k,:);
end
% Linear solver, Post-processing, classic ANN
[~,~,t_lin,~,~,~] = solveAndNormalize(0,imdl_lin,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_ann,~,~,~] = solveAndNormalize(0,imdl_ann,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_pp,~,~,~] = solveAndNormalize(0,imdl_pp,bk,do_parallel,[t_vh;t_vi]');
[~,~,t_nonl,~,~,~] = solveAndNormalize(0, imdl_nonl, bk, do_parallel,[t_vh;t_vi]');
% Free space
vh_test=t_vh; vi_test=t_vi; img_lin=t_lin; 
img_pp=t_pp; img_ann=t_ann; img_nonl=t_nonl;
clear t_vh t_vi t_lin t_pp t_ann t_nonl stim msel;
fmdl = rmfield(fmdl,f_to_rem); 
imdl_lin.fwd_model = fmdl; 
imdl_pp.fwd_model = fmdl;
imdl_ann.fwd_model = fmdl;
imdl_nonl.fwd_model = fmdl;
if distort && do_simulation; 
    fmdl_sd = arrayfun(@(x) rmfield(x,f_to_rem), fmdl_sd); 
end

% Less images, keep fmdl so we can easily display images
% [img_lin(:).fwd_model] = deal(fmdl); [img_pp(:).fwd_model] = deal(fmdl);
% [img_ann(:).fwd_model] = deal(fmdl); [img_nonl(:).fwd_model] = deal(fmdl);
if distort && do_simulation
	c = mat2cell(fmdl_sd,ones(n_imgs_test,1),1);
    [inh_t(:).fwd_model] = deal(c{:});
	inh_t = inh_t(1,:); % They are all the same
else
	[inh_t(:).fwd_model] = deal(fmdl);
end

%% 6 Display result
n_pts = 256;
h = fspecial('average',7);
% pos_targs = [area_test{1}(5,:); -area_test{1}(6,:)]*(n_pts/2)+(n_pts/2);
% sz_targs = area_test{1}([1 2],:)*(n_pts/2);
tit_fig = {'Initial image', 'One step GN','PDIPM','ANN','one-step GN + ANN'};
if ~do_simulation
    tit_fig = tit_fig(2:end);
end
	
for k= 1:1:n_imgs_test
    if n_targs == 1
        pos_targs = [area_test(k,5); -area_test(k,6)]*(n_pts/2)+(n_pts/2)';
        sz_targs = area_test(k,[1 2])*(n_pts/2)';
    else
        pos_targs = [area_test{k}(5,:); -area_test{k}(6,:)]*(n_pts/2)+(n_pts/2);
        sz_targs = area_test{k}([1 2],:)*(n_pts/2);
    end
	to_disp = [repmat(inh_t(1,k),1,1) img_lin(1,k) img_nonl(1,k) ...
											img_ann(1,k) img_pp(1,k)];
	[to_disp(:).fwd_model] = deal(fmdl);
	if ~do_simulation
        to_disp = to_disp(2:end);
    end
	[~]=show_multi_fem(to_disp,'abscissa',tit_fig);
	[to_disp(:).calc_colours] = deal(struct('npoints',n_pts));
	[to_disp(:).calc_slices] = deal(struct('filter',h));
	[~,sub_plots]=show_multi_fem(to_disp,'abscissa',tit_fig,...
							'slice',[Inf Inf 1]);
	for l=1:length(to_disp)
		axes(sub_plots(l)); hold on;
		drawEllipse([pos_targs(:,1); sz_targs(:,1)]','g','LineWidth',1.5);
		drawEllipse([pos_targs(:,2); sz_targs(:,2)]','g','LineWidth',1.5);
	end
	
end
							
%% 7 Errors
clear RMS* AR PE RES SD RNG dRES;
inh_err = inh_t; [inh_err(:).fwd_model] = deal(fmdl);
inh_t_n = normalize_mdl_dir(inh_err,'higher');
img_lin_n = normalize_mdl_dir(img_lin,'higher');
img_ann_n = normalize_mdl_dir(img_ann,'higher');
img_pp_n = normalize_mdl_dir(img_pp,'higher');
img_nonl_n = normalize_mdl_dir(img_nonl,'higher');
if ~do_simulation
    % Actual position in phantom
	xyzr = repmat( {[area_test{1}([5 6],:); [1 1]; area_test{1}(1,:)]}, 1 , n_imgs_test);
else
%     xyzr = [area_test(:,[4 5]), ones(n_imgs_test,1), area_test(:,1)]';
    for k = 1:1:n_imgs_test
        xyzr{k} = [area_test{k}([4 5],:); ones(1,n_targs); area_test{k}(1,:)];
    end
end
imgs_GREIT = [inh_err; img_lin; img_nonl; img_ann; img_pp];
for k = 1:1:n_imgs_test
	RMS_lin(:,k) = calc_errorRMS(inh_t_n(k), img_lin_n(:,k), 'area');
	RMS_ann(:,k) = calc_errorRMS(inh_t_n(k), img_ann_n(:,k), 'area');
	RMS_pp(:,k) = calc_errorRMS(inh_t_n(k), img_pp_n(:,k), 'area');
	RMS_nonl(:,k) = calc_errorRMS(inh_t_n(k), img_nonl_n(:,k), 'area');
	RMS_err = [RMS_lin; RMS_nonl; RMS_ann; RMS_pp];
	[~, AR{k}, PE{k}, dRES{k}, SD{k}, RNG{k}] = eval_GREIT_multi_targ( imgs_GREIT(2:end,k), xyzr{k});
	PE_sep(k,:,1) = PE{k}(1,:); PE_sep(k,:,2) = PE{k}(2,:);
	dRES_sep(k,:,1) = dRES{k}(1,:); dRES_sep(k,:,2) = dRES{k}(2,:);
	SD_sep(k,:,1) = SD{k}(1,:); SD_sep(k,:,2) = SD{k}(2,:);
end
% dRES = abs(RES(1,:) - RES(2:end,:));

disp('RMS Errors')
disp('One step GN   PDIPM      ANN      Post processing')
disp(RMS_err')
disp('')
disp('Position errors (For both targets)')
disp('One step GN   PDIPM      ANN      Post processing')
disp(PE_sep)
disp('')
disp('Difference of resolution (For both targets)')
disp('One step GN   PDIPM      ANN      Post processing')
disp(dRES_sep)
disp('')
disp('Shape deformation (For both targets)')
disp('One step GN   PDIPM      ANN      Post processing')
disp(SD_sep)

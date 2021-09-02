function data = fwd_solve_add_noise(fwd_model, img)
% FWD_SOLVE_ADD_NOISE: data= fwd_solve_add_noise( img)
% Forward solve the problem, and add noise
% Input:
%    img       = image struct
% Output:
%    data = measurements struct, with noise
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)
% Options for noise
%	img.fwd_model.add_noise.noise = noise model, determined from real
%	measurements (funciton est_SNR_multiple_rec)
%	img.add_noise.noise_db, SNR of noise, in dB, if noise model is not
%	specified
%	img.fwd_model.add_noise.solve = forward solver to call (before adding
%	noise)
%
%	(C) 2015/07/23 Sebastien Martin


% correct input parameters if function was called with only img
if nargin == 1
    img = fwd_model;
elseif strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
    warning('EIDORS:DeprecatedInterface', ...
        ['Calling FWD_SOLVE_ADD_NOISE with two arguments is deprecated and will cause' ...
        ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

n_elecs = length(fwd_model.electrode);
is_open = isOpenDomain(fwd_model);

if isfield(fwd_model,'parameters') && isfield(fwd_model.parameters,'isLungs')
    is_lungs = fwd_model.parameters.isLungs;
else
    is_lungs = false;
end

try
    noise_mdl = fwd_model.add_noise.noise;
catch
    try
        noise_mdl = get_default_noise(n_elecs,is_open,is_lungs);
    catch
        try
            noise_db = fwd_model.add_noise.noise_db;
        catch
            error('You should specify a noise model or at least a noise value, in dB');
        end
        noise_mdl = repmat(noise_db, size(fwd_model.meas_select));
    end
end
try
    filt = fwd_model.add_noise.filter; % Filter data
catch
    filt = true;
end
try
    fwd_sol = fwd_model.add_noise.solve; % Actual forward solver
catch
    fwd_sol = 'eidors_default';
end

img.fwd_model.solve = fwd_sol;
data = fwd_solve(img); % Solve the forward problem with the forward solver

% noisy = data.meas;
noisy = EIT_add_wgn_meas(data.meas, noise_mdl(fwd_model.meas_select,:), filt, n_elecs);

% Adjust the data structure to return
data.name= [data.name, ' with noise'];
data.meas = noisy;

end
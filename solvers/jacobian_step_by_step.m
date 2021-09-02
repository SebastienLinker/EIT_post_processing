function [ J ] = jacobian_step_by_step( fwd_model, img )
%JACOBIAN_STEP_BY_STEP Determines the Jacobian step by step, i.e. almost line
%	by line
%   Jacobian calculation may be slow because of large memory requirements
%	Uses a Divide And Conquer approach to overcome Matlab's weakness
%
%	(C) 2015/08/20 Sebastien Martin

if nargin == 1
    img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
    warning('EIDORS:DeprecatedInterface', ...
        ['Calling JACOBIAN_STEP_BY_STEP with two arguments is deprecated and will cause' ...
        ' an error in a future version. First argument ignored.']);
end

% Set correct fields
try img.fwd_model.jacobian = img.fwd_model.jacobian_SbS.jacobian_func;
catch; img.fwd_model.jacobian = 'eidors_default';
end

% Predict the size of the Jacobian Matrix
is2D = (size(img.fwd_model.nodes,2)==2);
n_meas = calc_n_meas(img.fwd_model);
% Size can also depends on the nodes, but we just want an estimate here
J_sz = n_meas * size(img.fwd_model.elems,1) * 8; % Assumes double precision
if (~is2D) && (J_sz>5e9) % 5GB
    issingle = true;
else
    issingle = false;
end


n_stims = size(img.fwd_model.stimulation,2);
if n_stims <= 150
    n_steps = 1;
else
    n_steps = n_stims/150; % How to divide the stimulations
end

% Actually gets the Jacobian
J = jacobian_SbS_internal(img, n_stims, n_steps, n_meas, issingle);

end

function [J] = jacobian_SbS_internal(img, n_stims, n_steps, n_meas, issingle)

sz_step = ceil(n_stims/n_steps);
n_MpStim = n_meas/n_stims; %measurements per stimulation
sz_J = [n_meas size(img.fwd_model.elems,1)];
if isfield(img.fwd_model,'coarse2fine')
    if sz_J(2)==size(img.fwd_model.coarse2fine,1); sz_J(2)=size(img.fwd_model.coarse2fine,2);
    else sz_J(2) = size(img.fwd_model.coarse2fine,1);
    end
end

k=1;
n_MpStep = sz_step*n_MpStim; % # measurements per step
while (k*sz_step < n_stims)
    tmp_img = img;
    tmp_img.fwd_model.stimulation = img.fwd_model.stimulation( ...
        ((k-1)*sz_step+1):1:(k*sz_step) );
    tmp_J = calc_jacobian(tmp_img);
    if k==1
        sz_J(2)=size(tmp_J,2);
        if issingle
            J = zeros(sz_J(1), sz_J(2), 'single');
        else
            J = zeros(sz_J(1), sz_J(2));
        end
    end
    if issingle; tmp_J = single(tmp_J); end
    J( ((k-1)*n_MpStep +1):1:(k*n_MpStep) ,:) = tmp_J;
    k=k+1;
end

if ((k-1)*sz_step < n_stims) %Last samples
    tmp_img = img;
    tmp_img.fwd_model.stimulation = img.fwd_model.stimulation( ...
        ((k-1)*sz_step +1):1:end );
    tmp_J = calc_jacobian(tmp_img);
    if k==1
        sz_J(2)=size(tmp_J,2);
        if issingle
            J = zeros(sz_J(1), sz_J(2), 'single');
        else
            J = zeros(sz_J(1), sz_J(2));
        end
    end
    if issingle; tmp_J = single(tmp_J); end
    J( ((k-1)*n_MpStep +1):1:end ,:) = tmp_J;
end

end

% Found in file calc_meas_icov.m
function n_meas = calc_n_meas( fwd_model )
n_meas = 0;
for i= 1:length(fwd_model.stimulation );
    n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern,1);
end
end
function img= inv_solve_diff_GN_one_step_no_RM( inv_model, data1, data2)
% INV_SOLVE_DIFF_GN_ONE_STEP inverse solver using approach of Adler&Guardo 1996
% img= inv_solve_diff_GN_one_step( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: inv_solve_diff_GN_one_step.m 3777 2013-03-20 15:30:44Z bgrychtol $

% EDIT: assume the reconstruction matrix has been calculated before, so we
% can save time

% Important: This is different from solve_use_matrix, in a sense that we
% guarantee the result to have double precision (in case the matrix
% contains single precision)

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

% Noise removal
try
    noise_th = inv_model.inv_solve.noise_th;
    dv = remove_noise(dv, noise_th);
catch
end

if isfield(inv_model,'RM')
    RM = inv_model.RM;
elseif isfield(inv_model,'parameters') && isfield(inv_model.parameters,'RM')
    RM = inv_model.parameters.RM;
elseif isfield(inv_model,'inv_solve') && isfield(inv_model.inv_solve,'RM')
    RM = inv_model.inv_solve.RM;
else
    warning('EIDORS:lateRM','Unknown Reconstruction Matrix. For better performance, RM should NOT be computed here');
    eidors_msg('Attempt late computation of RM',2);
    RM = get_RM_GN_one_step(inv_model);
end

sol = RM * dv;
sol = double(sol);

img = data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by inv_solve_diff_GN_one_step_no_RM';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
img = data_mapper(img,1);
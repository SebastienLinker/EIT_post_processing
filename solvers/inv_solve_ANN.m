function img= inv_solve_ANN( inv_model, varargin)
% INV_SOLVE_ANN inverse solver using ANN
% img= inv_solve_ANN( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2015/05/05 Sebastien Martin
%

if length(varargin)==2
    dv = calc_difference_data(varargin{1}, varargin{2}, inv_model.fwd_model);
else
    dv = varargin{1};
    try
        if strcmp(dv.type,'data')
            dv = dv.meas;
        end
    catch
    end
end

% Noise removal
try
    noise_th = inv_model.inv_solve.noise_th;
    dv = remove_noise(dv, noise_th);
catch
end

try
    ann = inv_model.parameters.ann;
catch
    ann = inv_model.inv_solve.ann;
end
sol = ann(dv);

img = data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by inv_solve_ANN';
try
    do_nodes = inv_model.(mfilename).do_nodes;
catch
    do_nodes = false;
end %Backward compatibility
if do_nodes
    img.node_data = sol;
    try
        n2e = inv_model.(mfilename).n2e;
    catch
        [~, n2e] = node_elem_mapper(inv_model);
    end
    img.elem_data = n2e * img.node_data;
else
    img.elem_data = sol;
end
img.fwd_model= inv_model.fwd_model;
img = data_mapper(img,1);

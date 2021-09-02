function [ imgs ] = smoothFEMdata( imgs, smooth_type, varargin )
%SMOOTHFEMDATA Apply smoothing filter on the conductivity distribution
%   The algorithm should be adapted to FEM

if isnumeric(smooth_type)
    varargin{1} = smooth_type;
    smooth_type = 'intuitive';
end
smooth_params = process_args(smooth_type,varargin{:});
n_imgs = numel(imgs);

% Make sure all the fmdls share the same structure
if ((length(imgs)>1) && (~isequal( imgs(:).fwd_model )))
    chk_elems = cell(1,n_imgs);
    for k=1:1:n_imgs
        chk_elems{k} = imgs(k).fwd_model.elems;
    end
    if ~isequal(chk_elems{:})
        error('All of the images should have the same model');
    end
end
fmdl = imgs(1).fwd_model;
[n2e, e2n] = node_elem_mapper(fmdl);
do_nodes = isfield(imgs(1),'node_data') && ~isfield(imgs(1),'elem_data');

if do_nodes
    img_data = [imgs(:).node_data];
    mat_smooth = n2e*e2n;
else %work on elements data
    img_data = [imgs(:).elem_data];
    if isfield(fmdl,'coarse2fine') && size(fmdl.coarse2fine,2)==size(img_data,1)
        img_data = fmdl.coarse2fine * img_data;
    end
    mat_smooth = e2n*n2e;
end

switch smooth_type
    case 'intuitive'
        for k = 1:1:smooth_params.n_it
            img_data = mat_smooth * img_data;
        end
    case 'average'
        adj_mat = build_adj_mat(fmdl, smooth_params.sz);
        for k = 1:1:n_imgs
            old_data = img_data(:,k);
            new_data = zeros(size(old_data));
            for l=1:1:length(fmdl.elems)
                new_data(l) = mean(old_data( adj_mat(l, adj_mat(l,:)~=0)));
            end
            img_data(:,k) = new_data;
        end
    otherwise
        error('This smoothing function is not implemented');
end

% img_data = mat2cell(img_data',ones(n_imgs,1),size(mat_smooth,1));
img_data = mat2cell(img_data, size(mat_smooth,1), ones(1,n_imgs));
if do_nodes
    [imgs(:).node_data] = img_data{:};
else
    [imgs(:).elem_data] = img_data{:};
end

end

function adj_mat = build_adj_mat(fmdl, sz)
pp = fwd_model_parameters(fmdl);
is2D = size(fmdl.elems,2)==3;
n_elems = size(fmdl.elems,1);
n_it = floor(sz/2 -1);
max_conn = sum(3*2.^(0:1:n_it));
adj_mat = zeros(n_elems,max_conn);
for k=1:1:n_elems
    tmp = find_adjoin(k, pp.ELEM, n_it);
    adj_mat(k,1:length(tmp)) = tmp;
end
end

% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM, n_it)
nn= ELEM(:,ee);
[d,e]= size(ELEM);
ss = false(size(ELEM));
for i=1:d
    ss= ss | ELEM==nn(i);
end
curr_elems= find(sum(ss,1)==d-1);
elems = curr_elems;
if n_it>=1
    for k=1:1:length(curr_elems)
        elems = [elems find_adjoin(curr_elems(k), ELEM, n_it-1)];
    end
    elems = unique(elems);
    elems( elems == (find(sum(ss,1)==d))) = [];
end
end

function params = process_args(smooth_type, varargin)
if strcmp(smooth_type,'intuitive')
    params.n_it = varargin{1};
elseif strcmp(smooth_type,'average')
    params.sz = varargin{1};
end
end
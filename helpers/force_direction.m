function [ imgs, switched ] = force_direction( imgs, dir, bk )
%SWITCH_DIRECTION Force direction of the target object
%   If the target has lower conductivity than background, switch to higher
%	If it has higher conductivity, switch to lower
%	Input
%		imgs: EIDORS images
%		dir: -1 (or 'lower') or 1 (or 'higher')
%	For display purpose only

debug = false;
n_imgs = length(imgs);
dir = process_input( n_imgs, dir );
switched = false(1,n_imgs);

for k=1:1:n_imgs
    img = imgs(k);
    % Estimate background and direction
    if nargin<3; bk = median(img.elem_data); end
    mini = min(img.elem_data); maxi = max(img.elem_data);
    if abs(mini-bk)>abs(maxi-bk); act_dir = -1; else; act_dir = 1; end
    % Verify and change direction if necessary
    if act_dir~=dir(k)
        if debug; figure; subplot(1,2,1); show_fem(img,1); title('Original'); end
        img = switch_dir(img, bk);
        if debug; subplot(1,2,2); show_fem(img,1); title('New direction'); end
        switched(k) = true;
    end
    imgs(k) = img;
end

end

function img = switch_dir(img, bk)
img.elem_data = -(img.elem_data-bk) + bk;
end

function dir = process_input( n_imgs, dir )
% Deal with char
if ischar(dir)
    if strcmp(dir,'lower') || strcmp(dir,'below')
        dir = -1;
    elseif strcmp(dir,'higher') || strcmp(dir,'above')
        dir = 1;
    else
        error('Unknown text');
    end
end
% Direction
if length(dir)==1; dir = repmat(dir,n_imgs,1); end
if length(dir)~=n_imgs; error('Unknown input format, variable dir'); end
end
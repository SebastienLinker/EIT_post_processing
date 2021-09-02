function [ normalized_set ] = normalize_mdl_dir( img_set, direction )
%NORMALIZE Normalize the data from 0 to 2
%   NOTE: This is different from mdl_normalize provided by EIDORS. This
%   function normalizes the conductivities, while EIDORS's function is
%   designed to normalize the measurements from electrodes
%	@input img_set: A set of images to normalize
%	@input direction: 'lower','higher' or 'auto', Tells whether the target
%	has a lower or higher conductivity than the background (default 'auto')
%		'2sides' Allows a normalization from two sides, in case of multiple
%		targets
%
%	This is an amelioration of normalize_mdl
%	Estimate the background automatically
%
%	(C) Sebastien Martin

debug = false;
img_set = squeeze(img_set);
n_iter = size(img_set);

if ~exist('direction','var'); direction = 'auto'; end
if strcmp(direction,'2sides'); do_2S = true; else do_2S=false; end

for k=1:1:n_iter(1)
    for l=1:1:n_iter(2)
        img = img_set(k,l);
        if debug; figure; show_fem(img,1); title('Initial image'); end
        % First, ensure we have more than 2 different conductivity values
        img_data = get_img_data(img);
        n_conduc = unique( img_data );
        if (size(n_conduc,1) == 1)
            error('Homogeneous image, WTF do you want to normalize?');
        elseif (size(n_conduc,1)==2)
            y = normalize_2values( img_data );
        else
            y = img_data;
            if strcmp(direction,'auto') % Automatic estimate of the direction
                white_bk = est_bk( mk_image(img.fwd_model, y));
                high=max(y);
                low=min(y);
                if (abs(high-white_bk)) > (abs(white_bk-low));
                    direction='lower';
                else direction='higher';
                end
            end
            switch direction
                case 'higher'
                    % Change direction
                    y = max(y)-y;
                    if debug; norm_debug(img.fwd_model,y,'Changed direction'); end
                case 'lower'
                    y = y-min(y);
                case '2sides'
                    white_bk = est_bk( mk_image(img.fwd_model, y));
                    bel=(y<white_bk); abo=(y>white_bk); is_bk=(y==white_bk);
            end
            if do_2S
                y = y-white_bk;
                y(bel) = y(bel)*0.5/(-min(y));
                y(abo) = y(abo)*0.5/(max(y));
                y(is_bk) = 0;
                y = y+0.5;
            else
                white_bk = est_bk( mk_image(img.fwd_model, y));
                y = y-white_bk;
                y = y./abs(min(y));
                if strcmp(direction,'lower')
                    y = 1+y;
                else
                    y=-y;
                end
                if debug; norm_debug(img.fwd_model,y,'Between 0 and 1'); end
            end
        end
        img.calc_colours.ref_level = 0; img.calc_colours.clim = 1;
        normalized = img;
        normalized.elem_data=y;
        if debug; figure; show_fem(normalized,1); title('Normalized image'); end
        normalized_set(k,l) = normalized;
    end
end

end

% 3 different possibilities
% Normalize from 0 to 2, from 0 to 1, or from 1 to 2
function y = normalize_2values( elem_data )
y = -1* ones(size(elem_data));
y( find(elem_data(:)==max(elem_data)) ) = 1;
y( find(elem_data(:)==min(elem_data)) ) = 0;
end

function [] = norm_debug( fem, elems, txttitle)
img.elem_data = elems;
img.fwd_model = fem;
img.type = 'image';
figure; show_fem(img,1);
if ~isempty(txttitle); title(txttitle); end
end
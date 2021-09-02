function [ normalized_set ] = normalize_mdl( img_set, bk, direction )
%NORMALIZE Normalize the data from 0 to 2
%   NOTE: This is different from mdl_normalize provided by EIDORS. This
%   function normalizes the conductivities, while EIDORS's function is
%   designed to normalize the measurements from electrodes
%	@input img_set: A set of images to normalize
%	@input bk: Conductivity of the backgorund (default 1)
%	@input direction: 'lower','higher' or 'auto', Tells whether the target
%	has a lower or higher conductivity than the background (default 'auto')


n_iter = numel(img_set);
debug = false;

if ~exist('bk','var'); bk = repmat(1,n_iter,1); end;
if ~exist('direction','var'); direction = 'auto'; end;
if length(bk)==1; bk=repmat(bk,n_iter,1); end;
if ischar(bk) && strcmp(bk,'graphics'); est_bk_graph = true; clear bk; else est_bk_graph = false; end

for k=1:1:n_iter
    img = img_set(k);
    % First, ensure we have more than 2 different conductivity values
    img_data = get_img_data(img);
    n_conduc = unique( img_data );
    if (size(n_conduc,1) == 1)
        error('Homogeneous image, WTF do you want to normalize?');
    elseif (size(n_conduc,1)==2)
        y = normalize_2values( img_data );
    else
        % The code below will return NaN in case we only have two different
        % conductivity values
        s_ed = img_data;
        s_ed(isnan(s_ed)) = [];
        s_ed = sort(s_ed);
        e = length(s_ed);
        if (e==0)
            error('Can''t display. All values NaN. Is raw data 0?')
        end
        y = img_data;
        if est_bk_graph
            n_taken = floor(.35*length(y));
            tmp = sort(y); bk(k) = mean(tmp((n_taken+1):(end-n_taken),:));
        end
        below_bkgnd = find(img_data<bk(k));
        div_below = (max(img_data(below_bkgnd))-s_ed(1));
        above_bkgnd = find(img_data>bk(k));
        div_above = (s_ed(e)-min(img_data(above_bkgnd)));
        switch direction
            case ('auto')
                abo=false;
                if isempty(div_above)
                    div = div_below; abo=false;
                elseif isempty(div_below)
                    div = div_above; abo=true;
                else
                    div = max(div_above, div_below);
                    if (div_above>div_below); abo=true; else abo=false; end;
                end
            case 'lower'
                div = div_below; abo=false;
            case 'higher'
                div = div_above; abo=true;
        end
        if abo
            y(below_bkgnd) = ( img_data(below_bkgnd)-max(img_data(below_bkgnd)) )/div;
        else y(below_bkgnd)=( img_data(below_bkgnd)-s_ed(1))/div;
        end
        if abo
            y(above_bkgnd) = (img_data(above_bkgnd)-min(img_data(above_bkgnd)) ) / div;
        else
            y(above_bkgnd) = 1+ (img_data(above_bkgnd)-min(img_data(above_bkgnd)))/div;
        end
    end
    normalized = img;
    normalized.elem_data=y;
    if debug; figure; show_fem(normalized,1); title('Normalized image'); end
    normalized_set(k) = normalized;
end

end

% 3 different possibilities
function y = normalize_2values( elem_data )
y = -1* ones(size(elem_data));
y( find(elem_data(:)==max(elem_data)) ) = 1;
y( find(elem_data(:)==min(elem_data)) ) = 0;
end
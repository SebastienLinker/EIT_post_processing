function [ greit, AR, PE, dRES, SD, RNG ] = eval_GREIT_multi_targ( images, xyzr, chk_dir )
%EVAL_GREIT_MULTITARG Evaluate GREIT errors for multi target images
%   Divide the images into two parts by tracing a line equidistant from the
%   two CoGs, and evaluate GREIT error on each target
%
%	@input images: A N*1 or 1*N array of N EIDORS images
%	@input xyzr: The xyzr (or xyr) coordinates of each target
%		A N*1 or 1*N cell array, containing the coordinates of each CoG
%		Each cell should contain a 3*M or 4*M array with coordinates for M
%		targets
%	@output greit: The greit errors
%
%	(C) 2015/04/07: Sebastien Martin

debug = false;
if nargin==2; chk_dir = true; end

n_imgs = length(images);
if ~iscell(xyzr); xyzr = {xyzr}; end
if length(xyzr)==1; xyzr = repmat(xyzr,n_imgs,1); end

for k=1:1:n_imgs
    img = images(k);
    is2D = (size(img.fwd_model.nodes,2)==2);
    n_elems = length(img.elem_data);
    c_xyzr = xyzr{k};
    bk = median(img.elem_data);
    n_targs = size(c_xyzr,2);
    if n_targs > 2;	error('This case is not implemented as yet'); end
    if size(c_xyzr,1)==3; c_xyzr = [c_xyzr([1 2],:);[0 0];c_xyzr(3,:)]; end;
    
    % Generate the equidistant line
    [sl, sh] = equidistantLine( reshape(c_xyzr(1:2,:),4,1) ); %[slope, shift]
    [img_sep, maps, c_xyzr] = separateTargets(img,sl,sh, c_xyzr);
    
    for l=1:1:n_targs
        c_img = rmfield( img_sep(l), 'elem_data');
        % Fill in the image with background
        c_img.fwd_model = reduceMesh( c_img.fwd_model, maps(l,:));
        c_img.elem_data = img_sep(l).elem_data;
        c_bk = median(c_img.elem_data);
        % If target < bk, then reverse
        haslower = (abs(c_bk-min(c_img.elem_data)) > abs(c_bk-max(c_img.elem_data)));
        if length(c_img.elem_data) > 0 && haslower && chk_dir
            c_img.elem_data = -(c_img.elem_data-c_bk)+c_bk;
        end
        c_img.elem_data = c_img.elem_data - c_bk;
        c_img2 = rmfield( img_sep(l), 'elem_data');
        c_img2.elem_data(maps(l,:)) = c_img.elem_data;
        c_img2.elem_data(~maps(l,:)) = 0;
        if debug; figure; show_fem(c_img2,1); title('1 target'); end
        % Estimate GREIT errors
        if is2D
            c_greit(l,:) = eval_GREIT_fig_merit(c_img2,c_xyzr(:,l));
        else
            c_greit(l,:) = eval_GREIT_3D(c_img2,c_xyzr(:,l));
        end
    end
    greit{k} = c_greit;
    AR(:,k) = c_greit(:,1); PE(:,k) = c_greit(:,2); dRES(:,k) = c_greit(:,3);
    SD(:,k) = c_greit(:,4); RNG(:,k) = c_greit(:,5);
end

end

function fmdl = reduceMesh(fmdl, toKeep)
fmdl.electrode = [];
fmdl.stimulation = [];
fmdl.meas_select = [];
fmdl.elems( ~toKeep,:) = [];
[fmdl.boundary, nodes] = find_boundary( fmdl.elems );
fmdl.name = 'Reduced model: Contains one target only. Use for error calculation only';
end
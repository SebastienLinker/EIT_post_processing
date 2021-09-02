function [ img_sep, maps, NCoGs ] = separateTargets( img, a, b, CoGs )
%SEPARATETARGETS Separate the image into two, according to the equidistant
%line between the CoG of the two targets
%   Call this function after calling equidistantLine

debug = false;
estCoGs = exist('CoGs','var');

fcn_abo = inline(['y>x*',num2str(a),'+',num2str(b)],'x','y','z');
abo = elem_select( img.fwd_model, fcn_abo);
maps(1,:) = (abo>0.5);
img_sep(1) = mk_image(img.fwd_model,1);
img_sep(1).elem_data = img.elem_data(maps(1,:));
fcn_bel = inline(['y<x*',num2str(a),'+',num2str(b)],'x','y','z');
bel = elem_select( img.fwd_model, fcn_bel);
maps(2,:) = (bel>0.5);
maps(2,:) = (abo<0.5); % bug fix
img_sep(2) = mk_image(img.fwd_model,1);
img_sep(2).elem_data = img.elem_data(maps(2,:));

if debug;
    img_disp = mk_image(img.fwd_model,0.5);
    img_disp.elem_data( maps(1,:) ) = 0;
    img_disp.elem_data( maps(2,:) ) = 1;
    figure; show_fem(img,1); title('Initial image, 2 targets');
    figure; show_fem(img_disp,1); title('Separated in 2 parts');
end

% Get the correct CoG for this subimage
if estCoGs
    for k=1:1:size(CoGs,2)
        isabo(k) = fcn_abo( CoGs(1,k), CoGs(2,k), CoGs(3,k));
        isbel(k) = fcn_bel( CoGs(1,k), CoGs(2,k), CoGs(3,k));
    end
    NCoGs = [CoGs(:,isabo), CoGs(:,~isabo)];
end

end
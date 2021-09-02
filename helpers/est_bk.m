function [ bk ] = est_bk( imgs )
%EST_BK Estimate background conductivity of images
% Estimates the background of the target image
% Considers the size of the elements
%
% (C) 2016/03/22 Sebastien Martin

n_imgs = numel(imgs);
bk = zeros(size(imgs));
for k=1:1:n_imgs
    img = imgs(k);
    if isfield('bk_ROI',img)
        selected{1} = inline(['x>',img.bk_ROI(1)],'x','y','z');
        selected{2} = inline(['x<',img.bk_ROI(2)],'x','y','z');
        selected{3} = inline(['y>',img.bk_ROI(3)],'x','y','z');
        selected{4} = inline(['y<',img.bk_ROI(4)],'x','y','z');
        selected = elem_select(img.fwd_model,selected)>0.5;
    else
        selected = true(size( get_img_data(img)));
    end
    % 	bk = median( img.elem_data(selected) ); % Doesn't consider the size of the elements, it
    % is a gross assumption for unbounded applications, avoid this simple line
    area_vec = abs(calc_elements_area(img.fwd_model));
    mean_vec = get_img_data(img) .* area_vec(:);
    bk(k) = mean(mean_vec(selected))/mean(area_vec(selected));
end
end
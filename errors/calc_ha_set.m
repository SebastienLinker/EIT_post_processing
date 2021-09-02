% Code from calc_hm_set adapted for 3D EIT
function hm_img = calc_ha_set(img,frac)
% hm_img= CALC_HA_SET(img)
% hm_img is sz*sz*sz*Nimg. It is 1 inside the Half Ampl Set
% frac is the fraction of maximum (0.5 or 0.25)
% hm_img expects conductive changes. Use calc_hm_set(-img,frac) for non-c

% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id: calc_hm_set.m 4809 2015-03-29 11:55:10Z bgrychtol-ipa $

sz = size(img);
[x,y,z]=meshgrid(linspace(-1,1,sz(1)),linspace(-1,1,sz(2)),linspace(-1,1,sz(3)));
map = x.^2+y.^2<1.1;

hm_img = logical(zeros(size(img)));
for i=1:size(img,4)
    imi = img(:,:,:,i); imi= imi(map);
    
    hmi= logical(zeros(sz));
    hmi(map) = imi >= (max(imi) * frac);
    hm_img(:,:,:,i) = hmi;
end
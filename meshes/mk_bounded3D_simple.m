function [ fmdl ] = mk_bounded3D_simple()
%MK_BOUNDED3D_SIMPLE Make the mesh similar to our physical phantom
%   Phantom has 4 layers of 8 electrodes
%   Big elements so computation is faster, less memory demanding
%	Add the two targets into the mesh itself
%	(C) 2021/08/14 Sebastien Martin

debug = false;
elec_per_layer = 8;
n_layers = 4;
elec_diam = 2;
thres_probe = 2;
h = 3;
bott_mid_cyl = 0; % Bottom of the medium cylinder
top_mid_cyl = 14; % Top of the medium cylinder
model_rad = 7;

th = linspace(0,2*pi,elec_per_layer+1)'; th(end) = [];
cs = [ model_rad*cos(th), model_rad*sin(th);];
elec_pos = [  cs, 0*th+(elec_diam/2+thres_probe), cs, 0*th];
idx=1.0;
for i=1:n_layers-1
    elec_pos = [elec_pos; cs, 0*th+(idx*h+elec_diam/2+thres_probe), cs, 0*th];
    idx=idx+1.0;
end

ref_top = '9.0'; ref_bot = '5.0'; % Just for illustration
obj1 = ['0,0,',ref_top,';0,0,',ref_bot,';6.1'];
contour = ['0,0,',num2str(top_mid_cyl),';0,0,0;6.4'];

shape_str = ['solid incyl = cylinder (0,0,',num2str(top_mid_cyl),'; 0,0,', ...
    num2str(bott_mid_cyl),'; ',num2str(model_rad),'); \n ', ...
    'solid pl1 =  plane(0,0,',num2str(bott_mid_cyl),';0,0,-1);\n ' ...
    'solid pl2 =  plane(0,0,',num2str(top_mid_cyl),'; 0,0,1);\n ', ...
    'solid obj1 = cylinder(',obj1,');\n ', ...
    'solid pl3 =  plane(0,0,',ref_bot,';0,0,-1);\n ', ...
    'solid pl4 =  plane(0,0,',ref_top,'; 0,0,1);\n ', ...
    'solid pl6 =  plane(0,0,3; 0,0,-1);\n ', ...
    'solid pl8 =  plane(0,0,',num2str(top_mid_cyl-3),'; 0,0,1);\n ', ...
    'solid contour = cylinder(',contour,');\n ', ...
    'solid obj_contour = pl6 and pl8 and contour;\n ', ...
    'solid targ_obj1 = pl3 and pl4 and obj1; tlo targ_obj1 -maxh=',num2str(3),';\n ',...
    'solid targ_contour = pl6 and pl8 and contour and not targ_obj1; tlo targ_contour -maxh=',num2str(8),';\n ',...
    'solid mainobj = (pl1 and pl2 and incyl) and not obj_contour; tlo mainobj;\n '];

elec_shape=[2, 2, 1];
for i=1:elec_per_layer*n_layers; elec_obj{i} = 'incyl'; end

if debug
    fmdl = ng_mk_gen_models(shape_str,elec_pos,elec_shape,elec_obj);
    figure; show_fem(fmdl); title('3D model');
else
    [~,fmdl] = evalc('ng_mk_gen_models(shape_str,elec_pos,elec_shape,elec_obj)');
end

fmdl = redo_elec_pos(fmdl,elec_pos,elec_shape,elec_obj, 7);
fmdl.nodes = fmdl.nodes ./ 7;

if debug; figure; show_fem(fmdl); end

end
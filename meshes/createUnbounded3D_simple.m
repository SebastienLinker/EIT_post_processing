function [ tank, s_tank, n_elements ] = createUnbounded3D_simple( n_lay, small_mdl )
%CREATEUNBOUNDED3DPHANTOM Creates the tank: unbounded 3D, without any
% target object
% For testing only, the output model is very coarse
%	Inputs:
%		n_lay: number of layers
%		small_mdl: To reduce the size, later
%
%	Outputs:
%		tank: A large tank that can be used to solve forward and inverse problems
%		s_tank: A part of this small tank, good for display. Note that
%			using this mesh to solve inverse problem may lead to large
%			distortions due to the close boundaries.

% (C) 2016/05/11 Sebastien Martin

debug = false;
t_locref = 2;

%% Parameters
probe_d=0;
n_elec = 8;
probe_rad = 1;
ref_noROI = [10 1];

[b_m_cyl, t_m_cyl, t_h_cyl, b_l_cyl, mod_r, th_pb, h, el_sz] = get_tank_param();

max_distance = [0.1 6]+probe_rad; % maximal distance, in mm (farther, the mesh has larger elements)
max_distance( max_distance>small_mdl ) = [];
max_distance(end+1) = small_mdl;
cog_ref = [5.5 0; 0 0; 0 0]; % COG of refined regions
elem_size = [10000 10000 10000 10000]; %maximal of elements, in mm^3
elem_size = elem_size(1:1:length(max_distance));

%% Creates tanks
th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
cs = [ probe_rad*cos(th), probe_rad*sin(th)-probe_d;];
elec_pos = [  cs, 0*th+el_sz/2+th_pb, cs, 0*th];
idx=1.0;
for i=1:n_lay-1
    elec_pos = [elec_pos; cs, 0*th+(idx*h+el_sz/2)+th_pb, cs, 0*th];
    idx=idx+1.0;
end

max_distance(end+1) = mod_r;

shape_str = ['solid incyl  = cylinder (0,0,',num2str(t_m_cyl),'; 0,0,',num2str(b_m_cyl),...
    '; ', num2str(probe_rad),'); \n ', ...
    'solid farcyl0 = cylinder (0,0,',num2str(t_m_cyl),'; 0,0,',num2str(b_m_cyl),'; ',...
    num2str(max_distance(1)),'); \n ',...
    'solid topcyl = cylinder (0,0,',num2str(t_h_cyl),'; 0,0,',num2str(t_m_cyl),'; ',...
    num2str(max_distance(end)),'); \n ', ...
    'solid bottomcyl = cylinder (0,0,',num2str(b_m_cyl),'; 0,0,',num2str(b_l_cyl),'; ',...
    num2str(max_distance(end)),'); \n '];

shape_str = [ shape_str ...
    'solid pl1    =  plane(0,0,',num2str(b_m_cyl),';0,0,-1);\n ' ...
    'solid pl2    =  plane(0,0,',num2str(t_m_cyl),'; 0,0,1);\n ' ...
    'solid pl3    =  plane(0,0,',num2str(t_m_cyl),';0,0,-1);\n ' ...
    'solid pl4    =  plane(0,0,',num2str(t_h_cyl),'; 0,0,1);\n ' ...
    'solid pl5    =  plane(0,0,',num2str(b_l_cyl),';0,0,-1);\n ' ...
    'solid pl6    =  plane(0,0,',num2str(b_m_cyl),'; 0,0,1);\n ' ...
    'solid probe = incyl and plane(0,0,',num2str(b_l_cyl),'; 0,0,-1); \n ' ];

shape_str = [ shape_str ...
    'solid farcyl',int2str(2),' = cylinder (',num2str(cog_ref(2,1)),',',num2str(cog_ref(2,2)),...
    ',',num2str(t_m_cyl),...
    ';',num2str(cog_ref(2,1)),',',num2str(cog_ref(2,2)),',',num2str(b_m_cyl),';',...
    num2str(max_distance(3)),');\n '];

shape_str = [ shape_str ...
    'solid farcyl',int2str(3),' = cylinder (',num2str(cog_ref(3,1)),',',num2str(cog_ref(3,2)),...
    ',',num2str(t_m_cyl),...
    ';',num2str(cog_ref(3,1)),',',num2str(cog_ref(3,2)),',',num2str(b_m_cyl),';',...
    num2str(max_distance(4)),');\n '];

shape_str = [shape_str ...
    'solid mainobj= pl1 and pl2 and farcyl2 and not probe; tlo mainobj;\n '];


k = length(max_distance)-1;
shape_str_large = [shape_str ...
    'solid largecyl',int2str(k),' = farcyl',int2str(k),' and pl1 and pl2 ', ...
    'and not mainobj; tlo largecyl',int2str(k),'; \n ',...
    'solid topcylobj= pl3 and pl4 and topcyl; tlo topcylobj;\n ', ...
    'solid bottomcylobj= pl5 and pl6 and bottomcyl; tlo bottomcylobj;\n ' ];


elec_shape=[el_sz, el_sz, 0.1];
for i=1:n_elec*n_lay; elec_obj{i} = 'probe'; end

if debug
    tank = ng_mk_gen_models( shape_str, elec_pos, elec_shape, elec_obj )
    figure; show_fem(tank); title('3D model');
    s_tank = ng_mk_gen_models( shape_str, elec_pos, elec_shape, elec_obj )
    figure; show_fem(s_tank); title('3D model, smaller');
else
    [~, tank] = evalc('ng_mk_gen_models( shape_str, elec_pos, elec_shape, elec_obj )');
    [~, s_tank] = evalc('ng_mk_gen_models( shape_str, elec_pos, elec_shape, elec_obj )');
end

% c2f mapping
n_elements = size(tank.elems,1);
c2f = mk_coarse_fine_mapping(s_tank,tank);
s_tank.coarse2fine = c2f;


end

% 	model_rad: Radius of model (default 70)
% 	thres_probe: distance (mm) between the bottom of the FEM and the 1st layer (default 0.5)
% 	h: distance between layers (center to center), in mm (default 1)
% 	elec_diam: diameter of one electrode (mm) (default 0.4)
function [b_m_cyl, t_m_cyl, t_h_cyl, b_l_cyl, mod_r, th_pb, h, elec_sz] = get_tank_param()
b_m_cyl = 0; t_m_cyl = 10/2.5;
t_h_cyl = 7.5; b_l_cyl = -5;
mod_r = 70;
th_pb = 0.5/2.5;
h = 2/2.5;
elec_sz = 1/2.5;
end
function [ tank ] = mk_unbounded2D_phantom( probe_rad, n_elec )
%MK_UNBOUNDED2D_PHANTOM Create the large tank for unbounded 2D problems
%   This tank has a very large number of nodes, so it'll be better to use
%   it for forward problems only
%	@ input
%		probe rad: radius of the probe
%		n_elec: number of electrodes, typically 8 or 16
%		elec_diam: diameter of electrodes

%% Parameters
probe_d=0;
elec_diam = 0.4; % diameter of one electrode (mm)
h = 1; % distance between layers (center to center), in mm
layer=1; %2D model

max_distance = [0.85 0.95 1 1.2 1.3 1.5 1.7 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.5 5 6]; %maximal distance, in mm
elem_size = [0.05 0.05 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.25 10]; %[10 100 1000]; %maximal of elements, in mm^3
elem_size = elem_size(max_distance > probe_rad);
max_distance = max_distance(max_distance > probe_rad);
model_rad = 70; % maximal position
thres_probe = 20; % distance (mm) between the top of the FEM and the 1st electrode's layer
H = h*(layer-1)+elec_diam; % Heigh of the electrode layers, from 1st to last layer
bott_mid_cyl = -0.1; %Bottom of the medium cylinder
top_mid_cyl = 0; %Top of the medium cylinder

%% Call to netgen
th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
cs = [ probe_rad*cos(th), probe_rad*sin(th)-probe_d;];
elec_pos = [  cs, 0*th-((top_mid_cyl-bott_mid_cyl)/2), cs, 0*th];
idx=1.0;
for i=1:layer-1
    elec_pos = [elec_pos; cs, 0*th-(idx*h+elec_diam/2+thres_probe), cs, 0*th];
    idx=idx+1.0;
end
max_distance(end+1) = model_rad;
shape_str = ['solid incyl  = cylinder (0,0,',num2str(top_mid_cyl),'; 0,0,',num2str(bott_mid_cyl),'; ',num2str(probe_rad),') -maxh=1.0; \n ', ...
    'solid farcyl0 = cylinder (0,0,',num2str(top_mid_cyl),'; 0,0,',num2str(bott_mid_cyl),'; ',num2str(max_distance(1)),') -maxh=1.0; \n ' ];
for k = 1:1:length(max_distance)-1
    shape_str = [ shape_str ...
        'solid farcyl',int2str(k),' = cylinder (0,0,',num2str(top_mid_cyl),';0,0,',num2str(bott_mid_cyl),'; ',num2str(max_distance(k+1)),') -maxh=',num2str(elem_size(k)),';\n '];
end;
shape_str = [ shape_str ...
    'solid pl1    =  plane(0,0,',num2str(bott_mid_cyl),';0,0,-1);\n ' ...
    'solid pl2    =  plane(0,0,',num2str(top_mid_cyl),'; 0,0,1);\n ' ...
    'solid probe = incyl and plane(0,0,',num2str(bott_mid_cyl),'; 0,0,-1); \n ' ];
for k = length(max_distance)-1:-1:1
    shape_str = [ shape_str ...
        'solid largecyl',int2str(k),' = farcyl',int2str(k),' and pl1 and pl2 and not farcyl',int2str(k-1),'; tlo largecyl',int2str(k),'; \n '];
end
shape_str = [shape_str ...
    'solid mainobj= pl1 and pl2 and farcyl0 and not probe; tlo mainobj;\n ' ];

elec_shape=[elec_diam, top_mid_cyl-bott_mid_cyl, 0.1];
for i=1:n_elec*layer; elec_obj{i} = 'probe'; end

fmdl3D = ng_mk_gen_models( shape_str, elec_pos, elec_shape, elec_obj );
tank = change_to_2D(fmdl3D,{0});

end
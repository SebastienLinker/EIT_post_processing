function img=inv_solve_diff_pdipm_fast( inv_model, data1, data2)
% INV_SOLVE_DIFF_PDIPM inverse solver for difference data using Primal/Dual interior point method
% img= ab_pdipm( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
%  inv_model.inv_solve_diff_pdipm.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.inv_solve_diff_pdipm.norm_image 1 or 2 (DEFAULT 2)
%  inv_model.inv_solve_diff_pdipm.beta     (default 1e-6)
%
% Parameters:
%  max_iters =  inv_model.parameters.max_iterations (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2008 Andrea Borsic. License: GPL version 2 or version 3
% $Id: inv_solve_diff_pdipm.m 3661 2012-11-20 21:18:01Z bgrychtol $


pp= process_parameters(inv_model);

fwd_model= inv_model.fwd_model;

d=calc_difference_data( data1, data2, fwd_model);

%  img_bkgnd = calc_jacobian_bkgnd( inv_model );
%  J = calc_jacobian( img_bkgnd );

alpha=calc_hyperparameter( inv_model );

% L=calc_R_prior( inv_model );
W= calc_meas_icov( inv_model );
% if pp.norm_data==1
%   W = sqrt(W); % sW is in units of volts
% end

[J, L, LTL, JTW] = get_J_L_LTL_JTW(inv_model);

if     pp.norm_data==2 && pp.norm_image==2
    x= pdipm_2_2( J,W,alpha*L,d, pp);
elseif pp.norm_data==2 && pp.norm_image==1
    x= pdipm_2_1( J,W,alpha*L,d, pp);
elseif pp.norm_data==1 && pp.norm_image==2
    x= pdipm_1_2( J,alpha*L,d, pp, LTL, JTW);
elseif pp.norm_data==1 && pp.norm_image==1
    x= pdipm_1_1( J,W,alpha*L,d, pp);
end

% create a data structure to return
img.name = 'inv_solve_diff_pdipm_faster';
img.elem_data = x;
img.fwd_model = fwd_model;
end

function s= pdipm_2_2( J,W,L,d, pp)
[s]= initial_values( J, L, pp);

R = L'*L;
ds= (J'*W*J + R)\(J'*W*(d - J*s) - R*s);
s= s + ds;
end

function m= pdipm_1_2( J,L,d, pp, LTL, JTW)
[m,x,~,sz]= initial_values( J, L, pp);
m = zeros(size(m,1),size(d,2));
x = zeros(size(x,1),size(d,2));

I_M = speye(sz.M, sz.M);
for loop = 1:pp.max_iter
    % Define variables
    f = J*m - d;             F= spdiag(f);
    X= spdiag(x);
    e = sqrt(f.^2 + pp.beta);E= spdiag(e);
    
    % Define derivatives
    %       dFc_dm = (I_M - X*inv(E)*F)*J;
    dFc_dm = (I_M - X/E*F)*J;
    dFc_dx = -E;
    dFf_dm = LTL;
    dFf_dx = JTW;
    
    numerator = -[dFc_dm, dFc_dx; dFf_dm, dFf_dx];
    Ex = E*x; f_Ex = f-Ex; clear Ex;
    JTWx = JTW*x;
    LTLm = LTL*m;
    denominator = [f_Ex; JTWx + LTLm;];
    dmdx = numerator \ denominator;
    clear numerator denominator;
    
    dm =             dmdx(      1:sz.N,:);
    dx = x_update(x, dmdx(sz.N+(1:sz.M),:));
    
    m= m + dm; x= x + dx;
    %       loop_display(i)
    pp = manage_beta(pp);
end
end

function m= pdipm_2_1( J,W,L,d, pp);
[m,jnk,y,sz]= initial_values( J, L, pp);

I_D = speye(sz.D, sz.D);
for loop = 1:pp.max_iter
    % Define variables
    g = L*m;                 G= spdiag(g);
    Y= spdiag(y);
    s = sqrt(g.^2 + pp.beta);S= spdiag(s);
    
    % Define derivatives
    dFf_dm = 2*J'*W*J;
    dFf_dy = L';
    dFc_dm = (I_D - Y*inv(S)*G)*L;
    dFc_dy = -S;
    
    dmdy = -[dFf_dm, dFf_dy; dFc_dm, dFc_dy] \ ...
        [ J'*W*(J*m-d) + L'*y; g-S*y ];
    
    dm =             dmdy(      1:sz.N );
    dy = x_update(y, dmdy(sz.N+(1:sz.D)));
    
    m= m + dm; y= y + dy;
    %       loop_display(i)
    % debug([mean(abs([m,dm])), mean(abs([y,dy]))]);
    pp = manage_beta(pp);
end
end

function m= pdipm_1_1( J,W,L,d, pp);
[m,x,y,sz]= initial_values( J, L, pp);

I_M = speye(sz.M,sz.M);
I_D = speye(sz.D,sz.D);
Z_N = sparse(sz.N,sz.N);
Z_DM= sparse(sz.D,sz.M);
for loop = 1:pp.max_iter
    % Define variables
    g = L*m;                 G= spdiag(g);
    r = sqrt(g.^2 + pp.beta);R= spdiag(r); % S in paper
    Y= spdiag(y);
    
    f = J*m - d;             F= spdiag(f);
    e = sqrt(f.^2 + pp.beta);E= spdiag(e);
    X= spdiag(x);
    
    % Define derivatives
    As1 = Z_N;
    As2 = (I_M - X*inv(E)*F) * J;
    As3 = (I_D - Y*inv(R)*G) * L;
    Ax1 = J'*W;
    Ax2 = -E;
    Ax3 = Z_DM;
    Ay1 = L';
    Ay2 = Z_DM';
    Ay3 = -R;
    B1  = J'*W*x + L'*y;
    B2  = f - E*x;
    B3  = g - R*y;
    
    DD = -[As1,Ax1,Ay1; ...
        As2,Ax2,Ay2; ...
        As3,Ax3,Ay3] \ [B1;B2;B3];
    
    dm = DD(1:sz.N);
    dx = x_update(x, DD(sz.N +        (1:sz.M)) );
    dy = x_update(y, DD(sz.N + sz.M + (1:sz.D)) );
    
    m= m + dm;
    x= x + dx;
    y= y + dy;
    %       loop_display(i)
    % debug([mean(abs([m,dm])), mean(abs([x,dx])), mean(abs([y,dy]))]);
    pp = manage_beta(pp);
end
end

% fix matlab's stupid verbose spdiags function
function sM = spdiag(V)
lV = length(V);
sM = spdiags(V,0,lV,lV);
end

function [s,x,y,sz]= initial_values( J, L, pp);
[sz.M,sz.N] = size(J); % M measurements, N parameters
[sz.D  ] = size(L,1); % E edges
y= zeros( sz.D, 1 ); % dual var - start with zeros
s= zeros( sz.N, 1 ); % solution - start with zeros
x= zeros( sz.M, 1 ); % dual var - start with zeros
end

% abs(x + dx) must be <= 1
function dx = x_update( x, dx)
% can't have zeros
dx(dx==0) = eps;
% space to limits in direction of dx
sx = sign(dx);
clr = sx - x;
% how much to multiply by to get to limits
fac = clr./dx;
% choose min amount to get to limits
dx = bsxfun(@times, dx, min(abs(fac)));
%  dx = dx*0.9;
end

function debug(vals)
%  disp(vals)
end

function pp = manage_beta(pp);
pp.beta = pp.beta * pp.beta_reduce;
if pp.beta < pp.beta_minimum;
    pp.beta = pp.beta_minimum;
end
end

function pp= process_parameters(imdl);
try    pp.max_iter = imdl.parameters.max_iterations;
catch;  pp.max_iter = 10;
end

try    pp.min_change = imdl.parameters.min_change;
catch;  pp.min_change = 0;
end

try    pp.beta = imdl.inv_solve_diff_pdipm.beta;
catch;  pp.beta = 1e-6;
end

pp.beta_reduce = 0.2;
pp.beta_minimum= 1e-16;

try    pp.norm_data = imdl.inv_solve_diff_pdipm.norm_data;
catch;  pp.norm_data = 2;
end

try    pp.norm_image = imdl.inv_solve_diff_pdipm.norm_image;
catch;  pp.norm_image = 2;
end
end

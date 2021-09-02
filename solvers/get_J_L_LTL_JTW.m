function [ J, L, LTL, JTW ] = get_J_L_LTL_JTW( imdl, W )
%GET_LTL_JTW Precalculating LTL and JTW in order to make inverse solver
%inv_solve_diff_pdipm faster
%   Use solver inv_solve_diff_pdipm_fast in addition to this script

if isfield(imdl,'RtR_prior') && ischar(imdl.RtR_prior)
    if strcmp(imdl.RtR_prior,'eidors_default')
        imdl.RtR_prior = eidors_default('get','calc_RtR_prior');
    end
    imdl.RtR_prior = str2func(imdl.RtR_prior);
end
if isfield(imdl,'R_prior') && ischar(imdl.R_prior)
    if strcmp(imdl.R_prior,'eidors_default')
        imdl.R_prior = eidors_default('get','calc_R_prior');
    end
    imdl.R_prior = str2func(imdl.R_prior);
end

% J
if isfield(imdl.inv_solve,'J')
    J = imdl.inv_solve.J;
else
    img_bkgnd = calc_jacobian_bkgnd( imdl );
    J = calc_jacobian( img_bkgnd );
end
% L
if isfield(imdl.inv_solve,'L')
    L = imdl.inv_solve.L;
else
    L=calc_R_prior( imdl );
end
% LTL
if isfield(imdl.inv_solve,'LTL')
    LTL = imdl.inv_solve.LTL;
else
    LTL = eidors_obj('get-cache',imdl,'LTL');
    if isempty(LTL)
        LT = L';
        LTL = LT*L;
    end
end
% JTW
if isfield(imdl.inv_solve,'JTW')
    JTW = imdl.inv_solve.JTW;
else
    JTW = eidors_obj('get-cache',imdl,'JTW');
    if isempty(JTW)
        JT = J';
        if nargin==1
            W= calc_meas_icov( imdl );
        end
        JTW = JT*W;
    end
end
end
function [ out ] = show_fem_dir( img, dir, transp, varargin )
%SHOW_FEM_DIR Displays EIDORS FEM with a target in a specific direction
%   No matter what is the target's conductivity (compared to background),
%   we always display in the same direction, usually 'higher'
%	varargin arguments are transmit to EIDORS function show_fem
%
%	(C) 2015/05/07 Sebastien Martin
%

% Check direction
if nargin == 1
    try dir = img.parameters.dir;
    catch
        eidors_msg('@@@ Unspecified direction, assumes higher',3);
        dir = 'higher';
    end
end
img = force_direction(img,dir);
% Actually displays FEM
if nargin < 3
    out = show_fem(img,varargin{:});
else
    out = show_transparent_fem(img,transp,varargin{:});
end

end
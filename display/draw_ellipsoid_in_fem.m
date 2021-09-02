function [ h_ell ] = draw_ellipsoid_in_fem( varargin )
%DRAW_ELLIPSOID_IN_FEM Draw an ellipsoid target inside the current axes
%handle
%   Make sure we don't get something too big that escapes the FE model
%
%	(C) 2015/11/17 Sebastien Martin

options = {'FaceColor', 'g', 'linestyle', 'none'};
is_closed = false;
if isAxisHandle(varargin{1})
    h = varargin{1};
    axes(h);
    varargin = varargin{2:end};
else
    h = gca;
end
k=1;
while k<=length(varargin)
    if ischar(varargin{k}) && strcmpi(varargin{k},'Closed')
        if k<nargin
            varargin(k:end) = varargin(k+1:end);
        else
            varargin = varargin(1:k-1);
        end
        is_closed = true;
    else
        k = k+1;
    end
end

[x,y,z] = drawEllipsoid(varargin{:});

x_lim = xlim;
y_lim = ylim;
z_lim = zlim;
if ~is_closed
    x(x<x_lim(1)) = x_lim(1)-0.01;
    x(x>x_lim(2)) = x_lim(2)+0.01;
    y(y<y_lim(1)) = y_lim(1)-0.01;
    y(y>y_lim(2)) = y_lim(2)+0.01;
    z(z<z_lim(1)) = z_lim(1)-0.01;
    z(z>z_lim(2)) = z_lim(2)+0.01;
else
    x(x<x_lim(1)) = x_lim(1)+0.01;
    x(x>x_lim(2)) = x_lim(2)-0.01;
    y(y<y_lim(1)) = y_lim(1)+0.01;
    y(y>y_lim(2)) = y_lim(2)-0.01;
    z(z<z_lim(1)) = z_lim(1)+0.01;
    z(z>z_lim(2)) = z_lim(2)-0.01;
end

c_bar = findobj(h, 'tag','Colorbar');
h_ell = surf(h, x, y, z, options{:});

end
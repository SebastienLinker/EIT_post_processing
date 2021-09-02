function [ h ] = drawTriangle( CoG, varargin )
%DRAWTRIANGLE draws a triangle on the existing plots
%   Detailed explanation goes here
%
%	(C) 2015/07/04 Sebastien Martin

if iscell(CoG)
    h = drawTriangle( CoG{:}, varargin{:});
    return
end
x = CoG; y = varargin{1}; sz = varargin{2}; alpha = varargin{3};
if nargin>=5; varargin = varargin(4:end); else varargin = {}; end

high = sqrt(3)/2*sz;
ang = 2*pi/3;
[theta, R] = cart2pol(2/3*high, 0);
% Define the position of the points
[pts(1,1), pts(1,2)] = pol2cart(theta+alpha, R);
[pts(2,1), pts(2,2)] = pol2cart(theta+ang+alpha, R);
[pts(3,1), pts(3,2)] = pol2cart(theta-ang+alpha, R);
pts = pts + repmat([x,y],3,1); %shift

% Actually plot the lines
h = drawEdge( [pts(1,:) pts(2,:); ...
    pts(2,:) pts(3,:); ...
    pts(3,:) pts(1,:); ], varargin{:});
end
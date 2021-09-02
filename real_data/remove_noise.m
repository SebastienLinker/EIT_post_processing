function [ dv ] = remove_noise( dv, varargin )
%REMOVE_NOISE Remove noise from difference data dv
%   Detailed explanation goes here

[method, args] = parse_varargin(varargin{:});

switch method
    case 'zero'
        dv( abs(dv)>args.th ) = 0;
end

end

function [method, args] = parse_varargin(varargin)
if length(varargin) == 1 % Default
    varargin{2} = varargin{1};
    varargin{1} = 'zero';
end
method = varargin{1};
% Parsing other arguments
switch method
    case 'zero'
        args.th = varargin{2};
end
end
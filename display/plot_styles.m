function [ h ] = plot_styles( Y, varargin )
%PLOT_STYLES Plot with different line styles, line colors, and markers
%   MATLAB plot function plots the different lines with different colours,
%   but the same linestyle and no marker, which gives a poor visibility
%   for published plots when printing in black and white.
%	This function modifies the linestyles and markers as well.
%
%	Input parameters are similar to MATLAB's function plot. Please refers
%	to plot function for further documentation
%	By default, it calls the 'plot' function, you can change it to other
%	function, for instance: plot_styles( @semilogx, ax, Y, ... )
%	NOTE: You may use varargin to define the color, linestyles, or markers,
%	but these properties will be overridden.
%
%	Copyright: Sebastien Martin 2015/07/13
%
%	See also PLOT

% Can call other plotting function, such as semilogx, plotyy, ...
if isa(Y,'function_handle')
    fcn = Y;
    Y = varargin{1};
    varargin = varargin(2:end); % Should have at least 2 arguments in this case
else
    fcn = @plot;
end

if ishandle(Y)
    ax = Y;
    Y = varargin{1};
    if length(varargin)>1
        varargin = varargin(2:end);
    else varargin = {};
    end
else ax = gca;
end

if ~isempty(varargin) && isnumeric(varargin{1})
    X = Y;
    Y = varargin{1};
    if length(varargin)>1
        varargin = varargin(2:end);
    else varargin = {};
    end
else
    X = 1:1:size(Y,1);
end
sz_legend = size(Y,2);

% set colours
colours = {'b';'g';'r';'k';'c';'m'}; n_colours = length(colours);
if n_colours<sz_legend; colours = repmat(colours, ceil(sz_legend/n_colours), 1); end
colours = colours(1:1:sz_legend);
% set markers
markers = {'o';'+';'*';'s';'d';'v';'>';'<';'^';'v';'h'}; n_markers = length(markers);
if n_markers<sz_legend; markers = repmat(markers, ceil(sz_legend/n_markers), 1); end
markers = markers(1:1:sz_legend);
% set linestyles
linestyle = {'-';'--';'-.';':'}; n_lines = length(linestyle);
if n_lines<sz_legend; linestyle = repmat(linestyle, ceil(sz_legend/n_lines), 1); end
linestyle = linestyle(1:1:sz_legend);

% Actually plot
h = fcn(ax,X,Y, varargin{:});

% Change markers, linestyle, and colors
mk_new = markers(end-length(h)+1:1:end);
line_new = linestyle(end-length(h)+1:1:end);
col_new = colours(end-length(h)+1:1:end);
cellfun( @(x,y,z,a) set(x,'Marker',y,'LineStyle',z,'Color',a), ...
    mat2cell(h,ones(size(mk_new))), mk_new, line_new, col_new);

end

function [ fig, h ] = show_multi_fem( imgs, varargin )
%SHOW_MULTI_FEM Display several EIDORS images on the same figure using
%subplot
%   Purpose: Makes comparison easily
%	Inputs:
%		imgs: A set of EIDORS images, should be a vector or matrix
%		Additional parameters should be 2 variables in the form {'name',value}
%		abscissa: abscissas legend
%		ordinate: ordinate legend
%		letters: boolean, adds the letters (a), (b), ... on subplots (default true)
%		axis: Specific axis range
%		colorbar: Adds the colorbar next to each indivual plot (default false)
%		reorder
%		FE_set: To set properties of the figure obtained from show_fem()
%		force_dir: Force direction of the target ('lower' or 'higher' than
%			the background
%		add_draw: Add an additional drawing over the FEM. Use a cell array
%		{drawing function, function arguments} to specify the properties of
%		the plot
%		suptitle: Adds a title above all plots
%		slice: Show a single slice instead of the FE model
%	Output: fig, the figure handle

%	(C) 2015/05/12 Sebastien Martin

debug = false;

imgs_ = squeeze(imgs);

[opts, hhprops] =  process_varargin(size(imgs),varargin{:});

if opts.reorder; imgs = reshape(imgs,opts.new_format); end
img_sz = size(imgs);
letter = 'a';
h = [];

if ~opts.disp_elec
    for k = 1:1:numel(imgs); imgs(k).fwd_model.electrode = []; end
end

if ~debug; fig = figure; end;
% Plot the FEMs
for k = 1:1:img_sz(1)
    for l = 1:1:img_sz(2)
        if debug
            h(end+1) = figure; fig = h;
        else
            h(end+1) = subtightplot(img_sz(1),img_sz(2), (k-1)*img_sz(2)+l); axis off;
        end
        if opts.slice.do_slice
            hh = show_slices(imgs(k,l),opts.slice.slice);
            % 			axis(opts.axis(l,:));
        elseif isempty(opts.force_dir)
            hh = show_transparent_fem(imgs(k,l),opts.transparency, opts.colorbar);
            axis(opts.axis(l,:));
        else
            hh = show_fem_dir(imgs(k,l),opts.force_dir,opts.transparency, opts.colorbar);
            axis(opts.axis(l,:));
        end
        if ~isempty(hhprops); set(hh,hhprops{:}); end
        
        if opts.add_draw
            hold on;
            opts.drawing.draw_fcn(opts.drawing.draw_params{:});
        end
        
        if opts.letters
            text(0.0,1.0,['(',letter,')'],'Units', 'Normalized', 'VerticalAlignment', 'Top', ...
                'FontName','Times New Roman','FontSize',10);
            if letter~='z';	letter = char(letter+1); else letter = 'A'; end
        end
    end
end

% Plot the text titles for asbcissas and ordinates
% Complies to IEEE standards
if ~debug
    for m=1:1:img_sz(2)
        % 		subtightplot(img_sz(1),img_sz(2), m); % 1st row
        axes( h(m) );
        % 		title( sprintf(opts.abscissa.txt{m})); %, 'FontName','Times New Roman','FontSize',10);
        text(0.5,1.1-0.03*opts.abscissa.n_lines, ...
            sprintf(opts.abscissa.txt{m}), 'Units', 'Normalized', 'VerticalAlignment','Middle',...
            'HorizontalAlignment','Center','FontSize',11,'FontWeight','bold');
    end
    for m=1:1:img_sz(1)
        % 		subtightplot(img_sz(1),img_sz(2), (m-1)*img_sz(2)+1); % 1st row
        axes( h((m-1)*img_sz(2)+1) );
        txt_ord = text(0.12,0.5, opts.ordinate.txt{m}, 'Units', 'Normalized', 'VerticalAlignment', 'Bottom',...
            'HorizontalAlignment','Center', 'Rotation',90); %,'FontName','Times New Roman','FontSize',10);
        set(txt_ord,'FontSize',11,'FontWeight','bold');
    end
end
if opts.do_suptitle
    suptitle(opts.suptitle);
end

end

function [opts, hhprops] =  process_varargin(img_sz,varargin)
% default
opts.abscissa = struct('txt','','n_lines',1);
opts.ordinate = struct('txt','','n_lines',1);
opts.letters = true;
opts.axis = 'auto';
opts.colorbar = 0;
opts.transparency = 1; % opaque
opts.reorder = false;
opts.new_format = []; % with reorder option, how to order the images
opts.force_dir = '';
opts.add_draw = false;
opts.do_suptitle = false;
hhprops = {}; %empty cell array
opts.slice.do_slice = false;
opts.disp_elec = true;
% loop
n_params = length(varargin);
k=1;
while (k<=n_params)
    opt = varargin{k};
    if ~ischar(opt); error(['Function ',mfilename,': Parameter not understood']); end
    switch opt
        case 'abscissa'
            opts.abscissa.txt = varargin{k+1}; k=k+2;
            opts.abscissa.n_lines = 1;
        case 'ordinate'
            opts.ordinate.txt = varargin{k+1}; k=k+2;
            opts.abscissa.n_lines = 1;
        case 'letters'
            opts.letters = varargin{k+1}; k=k+2;
        case 'axis'
            opts.axis = varargin{k+1}; k=k+2;
        case 'colorbar'
            opts.colorbar = varargin{k+1}; k=k+2;
        case 'transparency'
            opts.transparency = varargin{k+1}; k=k+2;
        case 'reorder'
            opts.reorder = true; opts.new_format = varargin{k+1}; k=k+2;
            if prod(opts.new_format)~=prod(img_sz)
                error('Inconsistent size of array');
            end
        case 'FE_set'
            hhprops = [hhprops varargin{k+1}]; k=k+2;
        case 'force_dir'
            opts.force_dir = varargin{k+1}; k=k+2;
        case 'add_draw'
            opts.add_draw = true;
            opts.drawing.draw_fcn = varargin{k+1}{1};
            opts.drawing.draw_params = varargin{k+1}{2};
            k=k+2;
        case 'suptitle'
            opts.do_suptitle = true; opts.suptitle = varargin{k+1}; k=k+2;
        case 'slice'
            opts.slice.do_slice = true; opts.slice.slice = varargin{k+1}; k=k+2;
        case 'no_electrode'
            opts.disp_elec = false; k=k+1;
    end
end

% Correction
if size(opts.axis,1)==1; opts.axis = repmat(opts.axis,img_sz(2),1); end
opts.abscissa = chkAxisTitles(opts.abscissa,img_sz);
opts.ordinate = chkAxisTitles(opts.ordinate,img_sz);
end

function titles = chkAxisTitles(titles,n_imgs)
titles.txt = param2str(titles.txt);
% Use n_imgs(1 or 2) instead
if length(titles.txt)==1
    titles.txt = repmat(titles.txt,max(n_imgs),1);
end
% Check number of lines
titles.n_lines = max( cellfun(@(x) length(x), strfind(titles.txt,'\n') ))+1;
end
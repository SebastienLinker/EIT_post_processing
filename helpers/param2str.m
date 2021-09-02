function [ titles ] = param2str( titles )
%PARAM2STR Convert the EIDORS parameters into a string equivalent
%	Parameters can be a list of prior functions
%
% (C) 2015/07/08 Sebastien Martin

if ischar(titles); titles = {titles}; end
if iscell(titles)
    titles = cellfun(@(x) toString(x), titles,'UniformOutput',false);
elseif isvector(titles)
    titles = arrayfun(@(x) toString(x), titles, 'UniformOutput',false);
end

end

function str = toString(in)
if ischar(in); str = in; end
if isa(in,'function_handle')
    str = func2str(in);
    str = regexprep(str,'^prior_',''); % In case it is a prior
    str = strrep(str,'_',' ');
end
if isnumeric(in)
    str = num2str(in);
    str = regexprep(str,'^2$','ST-EIT'); % No CS in this case
    str = regexprep(str,'^0.','\\Delta\\alpha=0.');
end
end
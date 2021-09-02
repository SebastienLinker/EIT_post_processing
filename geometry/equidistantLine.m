function [ a, b ] = equidistantLine( CoGs )
%EQUIDISTANTLINE Return the equation of the line equidistant to two points
%   Detailed explanation goes here

is2D = (length(CoGs)==4);
if is2D
    pt1 = CoGs([1 2]);
    pt2 = CoGs([3 4]);
else
    pt1 = CoGs([1 2 3]);
    pt2 = CoGs([4 5 6]);
end

m_p = (pt1+pt2)/2; % Midpoint
% Slope of line that crosses the two points
if is2D
    slope_AB = (pt2(2)-pt1(2)) / (pt2(1)-pt1(1));
else
    error('Write code (o)(o)');
end
a = -1/slope_AB; % Slope of the equidistant line
% Equation y=ax+b
% m_p(2) = a*m_p(1)+b
b = m_p(2) - a*m_p(1);

end
function [xnorm,ynorm,min_range,max_range] = normalise(x,y)

% Normalise velocity vector in the range between -1 and 1. Original sign of
% the velocities is maintained.

% if ~(length(x) == length(y))
%     error('Input x,y must be of equal size')
% end

% merge zonal, meridional velocity
data = [x y];

% normalise while preserving sign
if abs(min(data)) > max(data)
    max_range = abs(min(data));
    min_range = min(data);
else
    max_range = max(data);
    min_range = -max(data);
end
norm_value = 2 .* data ./ (max_range - min_range);

xnorm = norm_value(1:length(x));
ynorm = norm_value(length(x)+1:end);

return